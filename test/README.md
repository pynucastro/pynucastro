# pyreaclib Unit Tests

The ```test``` directory contains unit testing for pyreaclib using nose. Run a unit test by executing the command:

```
$ nosetests -s
```

The ```test/standard``` subdirectory contains subdirectories for each problem in ```pyreaclib/examples``` with the standard output we expect pyreaclib to produce for that problem.

The ```runs``` directory is an empty directory (contents will be ignored by git) which holds the output of each unit test in a directory named according to the following pattern:

```
[year]-[month]-[day]--[hh]-[mm]-[ss]-[microseconds]
```

Within each such directory (eg. XYZ), ```test_reaclib.py``` will mirror the directory structure in ```test/standard``` while copying in the contents of the corresponding directories in ```pyreaclib/examples```. Next, ```test_reaclib.py``` runs the reaclib script in each of the directories in ```test/runs/XYZ``` and calls ```diff``` to compare each corresponding file in ```test/standard``` with the files produced in ```test/runs/XYZ```.

```test_reaclib.py``` expects the script file used to run the example to be named as follows for an example named ```ABC``` or ```ABC_DEF```: ```abc.py```

  

