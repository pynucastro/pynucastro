#!/usr/bin/env python

# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.
#
# M. Zingale (2012-03-21)

from __future__ import print_function
import sys
import re
import string
import os
import argparse
import subprocess

class Search(object):
    def __init__(self, prefix=None, program=None, vpath=None, files=None, output='dep'):
        self.prefix = prefix
        self.program = program
        self.vpath = vpath
        self.files = files
        self.output = output

        try:
            self.outfile = open(self.output, 'w')
        except:
            raise

        try:
            self.logfile = open(self.output + '.log', 'w')
        except:
            raise
        
        self.vpath_sources_f90 = []
        
        # keyed by module name, entry is full file path
        self.module_files = {}
        
        # keyed by source file base name, entry is list of dependency modules
        self.f90sources_depmods = {}

        # keyed by source file base name, entry is list of dependency source files (base names)
        self.f90sources_f90mods = {}

        # regular expression for ' use modulename, only: stuff, other stuff'
        # see (txt2re.com)
        self.use_re = re.compile("( *)(use)(\s+)((?:[a-z_][a-z_0-9]+))", 
                            re.IGNORECASE|re.DOTALL)
        self.module_re = re.compile("( *)(module)(\s+)((?:[a-z][a-z_0-9]+))",
                               re.IGNORECASE|re.DOTALL)
        self.module_proc_re = re.compile("( *)(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                                    re.IGNORECASE|re.DOTALL)
#        self.program_re = re.compile("( *)(program)(\s+)((?:[a-z_][a-z_0-9]+))", 
#                                re.IGNORECASE|re.DOTALL)

        self.vpath_to_f90()
        self.find_all_use_deps(self.files)
        self.modules_to_files()
        self.print_source_dependencies()
        self.print_sources()

    def lprint(self,s=''):
        """
        Write s to self.logfile, appending a newline.
        """
        self.logfile.write(str(s) + '\n')

    def oprint(self,s=''):
        """
        Write s to self.outfile, appending a newline.
        """
        self.outfile.write(str(s) + '\n')

    def vpath_to_f90(self):
        """
        Find all .f90 or .F90 source files in vpath, 
        adding their full paths to self.vpath_sources_f90
        """
        self.lprint('Run vpath_to_f90 ----------')
        for path in self.vpath:
            if os.path.isdir(path):
                self.lprint('* Looking in directory {}'.format(path))
                contents = os.listdir(path)
                self.lprint('* * Found contents:')
                self.lprint(contents)
                # Find all the .f90 source files
                for content in contents:
                    self.lprint('Checking {}'.format(content))
                    if os.path.isfile(os.path.join(path,content)):
                        base, ext = os.path.splitext(content)
                        if ext=='.f90' or ext=='.F90':
                            self.vpath_sources_f90.append(os.path.join(path, content))
                            self.lprint('* * Adding file {} to vpath_sources_f90'.format(content))
        self.lprint('End vpath_to_f90 ----------')

    def find_all_use_deps(self, files):
        """
        Find all the use dependencies in self.files, 
        filling self.f90sources_depmods
        """
        for file in files:
            self.find_use_dependencies(file)

    def find_use_dependencies(self, file):
        """
        Get the list of modules this file depends on and add 
        to the dictionary self.f90sources_depmods.

        Then call find_module_file on each of the dependency modules.
        """
        self.lprint('Run find_use_dependencies ----------')
        source_mod_dep = [] # List of modules this file uses
        f = open(file, "r")
        for line in f:
            # strip off the comments
            idx = string.find(line, "!")
            line = line[:idx]
            rebreak = self.use_re.search(line)
            if rebreak:
                source_mod_dep.append(rebreak.group(4)) # Module name
        self.lprint('* File {} depends on modules:'.format(file))
        self.lprint(source_mod_dep)
        self.f90sources_depmods[os.path.basename(file)] = source_mod_dep
        for module in source_mod_dep:
            self.find_module_file(module)
        self.lprint('End find_use_dependencies ----------')

    def modules_to_files(self):
        """
        Take the dependency modules and find dependency source files
        containing those modules.

        Translates self.f90sources_depmods to self.f90sources_f90mods
        """
        self.lprint('Run modules_to_files ----------')
        self.lprint('f90sources_depmods:')
        self.lprint(self.f90sources_depmods)
        for f90source in self.f90sources_depmods.keys():
            modulelist = self.f90sources_depmods[f90source]
            self.lprint('* File {} depends on:'.format(f90source))
            self.lprint(modulelist)
            filelist = []
            for module in modulelist:
                if module in self.module_files.keys():
                    filelist.append(os.path.basename(self.module_files[module]))
                else:
                    self.lprint('No file entry containing module {} in self.module_files!\n'.format(module))
            self.lprint(filelist)
            self.f90sources_f90mods[f90source] = filelist
        self.lprint('End modules_to_files ----------')

    def module_in_file(self, module, file):
        """
        Given a module and a source file, returns:
        * True, if module declaration is present in the file
        * False, otherwise
        """
        try:
            f = open(file, 'r')
        except:
            raise
        found = False
        for line in f:
            remod = self.module_re.search(line)
            procmod = self.module_proc_re.search(line)
            if remod and not procmod:
                found_module_name = remod.group(4)
                if module==found_module_name:
                    found = True
                    break
        f.close()
        return found

    def find_module_file(self, module):
        """
        Given a module, will search self.vpath_sources_f90 for it.
        When a file is found, adds it to self.module_files and stops looking.
        
        NOTE: This means it searches source files in the same order as
        the directories are listed in vpath!
        """
        if module in self.module_files.keys():
            self.lprint('Module {} already located in {}'.format(
                module,
                self.module_files[module]))
            return
        else:
            self.lprint('Run find_module_file ----------')
            self.lprint('* Looking for module {}'.format(module))
            for file in self.vpath_sources_f90:
                if self.module_in_file(module, file):
                    self.lprint('* Found module {} in file {}'.format(
                        module, file))
                    self.module_files[module] = file
                    # Add file to f90sources_depmods if not already there
                    if not file in self.f90sources_depmods.keys():
                        self.lprint('* File {} not in f90sources_depmods'.format(file))
                        self.find_use_dependencies(file)
                    self.lprint('End find_module_file ----------')
                    return
            # module not found
            self.lprint('Error, could not find module {} in vpath:\n'.format(module))
            self.lprint(self.vpath)
            self.lprint('End find_module_file ----------')

    def print_source_dependencies(self):
        """
        Print the file-file dependencies in self.f90sources_f90mods
        """
        self.lprint('Run print_source_dependencies ----------')
        for source in self.f90sources_f90mods.keys():
            sourcedeps = self.f90sources_f90mods[source]
            for sdep in sourcedeps:
                self.oprint(self.prefix + source.replace(".f90", ".o") +
                            ' : ' +
                            self.prefix + sdep.replace(".f90",".o"))
            self.oprint()
        self.oprint()
        self.lprint('End print_source_dependencies ----------')

    def print_sources(self):
        """
        Print the final list of f90sources
        """
        self.lprint('Run print_sources ----------')
        self.oprint("f90sources := {}".format(" ".join(
            self.f90sources_f90mods.keys())))
        self.lprint('End print_sources ----------')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--program",
                        help="name of program for which to prepare compilation dependencies",
                        default="")
    parser.add_argument("--vpath",
                        type=str,
                        default=None,
                        nargs="*",
                        help="This option takes an argument list consisting of the VPATH directories to search for source files.")
    parser.add_argument("--files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    parser.add_argument("--output", type=str, default='dep',
                        help='Name of the output file to be written. A separate log file will also be produced with ".log" appended to the output file name.')
    
    args = parser.parse_args()
    s = Search(args.prefix, args.program, args.vpath, args.files, args.output)




