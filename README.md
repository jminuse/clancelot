# clancelot
A set of computational chemistry Python libraries and tools developed by the Clancy Group, Cornell University.

## Installing
If you would like to install, simply clone the git repo, edit the install.py accordingly and run install.py

		`git clone https://github.com/jminuse/clancelot.git`

## Contributing
If added to the project and you want to contribute, you can do the following:

1. Request permission from James Stevenson
2. Setup an ssh key as described [here](https://help.github.com/articles/generating-ssh-keys/)
3. Decide where you would like to have the repo, change directories to it, and run:

		`git clone git@github.com:jminuse/clancelot.git`
4. Run the following to make sure you have your ssh id's stored:

		`ssh-add ~/.ssh/id_rsa`

		`ssh-add`
5. Whenever you want to add or delete a file (remember that you cannot add an empty folder) run:

	To Add:

		`git add PATH/TO/FILE/NEW_FILE_NAME`
	To Delete:

		`git rm PATH/TO/FILE/FILE_TO_DELETE`
6. Whenevery you want to add updated files to git, run:

		`git add -u`
7. When you are ready to add to the main repo, run:

		`git commit -m "Whatever message describing what you did"`
		`git push`
8. Finally, you may want to check if there are any changes you made prior to adding all updated files:
	The following checks which files are different to what is in your commit:
	
		`git status`
	The following will look at the last two changes of a file:
	
		`git blame PATH/TO/FILE/FILE_NAME`


TODO

1. Pull out NBS  
2. Subclass Exception to ClancelotException  
3. Sublclass DeprecationWarning to ClancelotDeprecationWarning  
4. Warnings & Exceptions & Deprecation file  
5. Documentation
  -
  -
  -
6. Switch to pep8 style (or some style)  
7. Quantum Espresso Support  
8. LAMMPs Support - abstract out input file  
9. Orca default routes? Best practices  
10. How to add new tests (commit hooks)  
11. 