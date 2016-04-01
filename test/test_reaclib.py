# Test the examples against the outputs in their corresponding ./standard/... subdirectories
from __future__ import print_function
import os
import sys
import importlib
import shutil
import datetime
import time
import difflib
import nose.tools

# Setup some directory names
dir_test = os.path.dirname(os.path.realpath(__file__))
dir_reaclib = os.path.split(dir_test)[0]
dir_examples = os.path.join(dir_reaclib, 'examples')
dir_runs = os.path.join(dir_test, 'runs')
dir_standard = os.path.join(dir_test, 'standard')

# Make a list of test cases from the directories in pyreaclib/test/standard
test_cases = []
for entry in os.listdir(dir_standard):
    if os.path.isdir(os.path.join(dir_standard, entry)):
        test_cases.append(entry)

# Use the current date and time to create a unique testrun directory name
dtnow = datetime.datetime.now()
sdtnow = dtnow.strftime("%Y-%m-%d--%H-%M-%S-%f")
dir_now = os.path.join(dir_runs, sdtnow)
# Make a directory for these tests in runs
os.mkdir(dir_now)
print('Using test suite directory {}'.format(dir_now))

def test_all():
    for tc in test_cases:
        print('\nTesting {}:\n'.format(tc))
        
        # Make dirs and copy script file
        dir_tc = os.path.join(dir_now, tc)
        dir_std_tc = os.path.join(dir_standard, tc)
        os.mkdir(dir_tc)
        example_path = os.path.join(dir_examples, tc)
        script_module_name = tc.split('_')[0].lower()
        for file in os.listdir(example_path):
            shutil.copy(os.path.join(example_path, file), dir_tc)
        os.chdir(dir_tc)

        # Add the test directory to path so I can import its script
        sys.path.insert(0, dir_tc)

        # Run test case
        mk = sys.modules.keys()
        mk_before = mk[:]
        tcmod = importlib.import_module(script_module_name)
        mk = sys.modules.keys()
        mk_after = mk[:]

        # Necessary if muliple test case paths use the same script name:
        ## Remove test case module and modules it imported from module list
        mk_to_del = list(set(mk_after)-set(mk_before))
        for k in mk_to_del:
            sys.modules.pop(k)

        ## Remove test directory from path
        sys.path.pop(0)
        
        # Compare with standard
        for file in os.listdir(dir_std_tc):
            f_tc = os.path.join(dir_tc, file)
            f_std_tc = os.path.join(dir_std_tc, file)
            try:
                lines_tc = open(f_tc, 'U').readlines()
                date_tc  = time.ctime(os.stat(f_tc).st_mtime)
            except:
                raise
            try:
                lines_std_tc = open(f_std_tc, 'U').readlines()
                date_std_tc = time.ctime(os.stat(f_std_tc).st_mtime)
            except:
                raise
            diff = difflib.unified_diff(lines_tc, lines_std_tc,
                                        f_tc, f_std_tc,
                                        date_tc, date_std_tc,
                                        n=3)
            diffreport = ''.join(diff)
            if diffreport.strip() != b'':
                print('')
                print(diffreport)
            yield nose.tools.assert_equals, diffreport.strip(), b''
    

