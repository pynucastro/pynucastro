#!/usr/bin/env python
# Given a date-time string XYZ corresponding to a directory in pyreaclib/test/runs,
# Replace the pyreaclib/test/standard files with the corresponding files in
# pyreaclib/test/runs/XYZ

from __future__ import print_function
import os
import shutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('datetime', type=str, nargs=1,
                    help='Date-time string corresponding to a directory in pyreaclib/test/runs. The files in pyreaclib/test/standard will be replaced by their corresponding files in pyreaclib/test/runs/datetime.')
args = parser.parse_args()

# Setup some directory names
dir_test = os.path.dirname(os.path.realpath(__file__))
dir_runs = os.path.join(dir_test, 'runs')
dir_standard = os.path.join(dir_test, 'standard')

# Make a list of test cases from the directories in pyreaclib/test/standard
test_cases = []
for entry in os.listdir(dir_standard):
    if os.path.isdir(os.path.join(dir_standard, entry)):
        test_cases.append(entry)

dir_run_replace = os.path.join(dir_runs, args.datetime[0])
if not os.path.isdir(dir_run_replace):
    print('ERROR: Supply a date-time string corresponding to a directory name within pyreaclib/test/runs')
    exit()
else:
    print('Using test suite directory {}'.format(dir_run_replace))

for tc in test_cases:
    print('\n----------------------------------------')
    print('Test case {}:'.format(tc))
    print('----------------------------------------\n')

    dir_tc = os.path.join(dir_run_replace, tc)
    dir_std_tc = os.path.join(dir_standard, tc)
    for file in os.listdir(dir_std_tc):
        print('* Replacing file ')
        print('* * {} '.format(os.path.join(dir_std_tc, file)))
        print('with')
        print('* * {}'.format(os.path.join(dir_tc, file)))
        print('')
        shutil.copy(os.path.join(dir_tc, file), dir_std_tc)

    

