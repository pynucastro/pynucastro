#!/usr/bin/env python

# automatically generate Makefile dependencies for Fortran 90 source.
#
# this will output all the dependency pairs amongst the source files.
#
# M. Zingale (2012-03-21)

import sys
import re
import string
import os
import argparse

def doit(prefix, program, files):

    # regular expression for ' use modulename, only: stuff, other stuff'
    # see (txt2re.com)
    use_re = re.compile("( *)(use)(\s+)((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    module_re = re.compile("( *)(module)(\s+)((?:[a-z][a-z_0-9]+))",
                           re.IGNORECASE|re.DOTALL)

    module_proc_re = re.compile("( *)(module)(\s+)(procedure)(\s+)((?:[a-z][a-z_0-9]+))",
                                re.IGNORECASE|re.DOTALL)
    
    program_re = re.compile("( *)(program)(\s+)((?:[a-z_][a-z_0-9]+))", 
                        re.IGNORECASE|re.DOTALL)

    # first parse the files and find all the module statements.  Keep a
    # dictionary of 'module name':filename.
    modulefiles = {}
    excludefiles = []
    for file in files:

        f = open(file, "r")
        
        for line in f:

            # strip off the comments
            idx = string.find(line, "!")
            line = line[:idx]

            rebreak = module_re.search(line)
            rebreak2 = module_proc_re.search(line)
            if rebreak and not rebreak2:
                modulefiles[rebreak.group(4)] = file
            reprogram = program_re.search(line)
            if reprogram and reprogram.group(4)!=program:
                excludefiles.append(file)
        f.close()

    # go back through the files now and look for the use statements.
    # Assume only one use statement per line.  Ignore any only clauses.
    # Build a list of dependencies for the current file and output it.
    f90sources = []
    for file in files:
        if file in excludefiles:
            continue

        f90sources.append(os.path.basename(file))
        f = open(file, "r")

        for line in f:

            # strip off the comments
            idx = string.find(line, "!")
            line = line[:idx]

            rebreak = use_re.search(line)
            if rebreak:
                print prefix+os.path.basename(file).replace(".f90", ".o"), ':', \
                    prefix+os.path.basename(modulefiles[rebreak.group(4)]).replace(".f90", ".o")

        f.close()
        print " "

    print "F90SOURCES := {}".format(" ".join(f90sources))
        
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix",
                        help="prefix to prepend to each dependency pair, e.g., for a build directory",
                        default="")
    parser.add_argument("--program",
                        help="name of program for which to prepare compilation dependencies",
                        default="")
    parser.add_argument("files", metavar="source files", type=str, nargs="*",
                        help="F90 source files to find dependencies amongst")

    args = parser.parse_args()
    doit(args.prefix, args.program, args.files)



