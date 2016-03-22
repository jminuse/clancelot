#!/bin/sh
#
# A hook script to test the code before it is committed.
# To enable this hook, copy this file to ".git/hooks/pre-commit" 
# and enable execution with chmod +x. This is done automatically
# by install.py. 
# This test does not stop the commit; it simply warns you if the 
# commit fails a test. To stop the commit, make this script 
# return a non-zero value. 

python unit_tests/pre_commit_tests.py

