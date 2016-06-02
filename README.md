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

1. Things to deprecate  
  [ ] In Utils:
  	[ ] frange, quat_to_mat, matvec, matmat
  	[ ] utils. elements_by_atomic_number  
  	[ ] files.inp_to_xyz  
  [ ] awk deprecated due to having in linter  
2. Things to adjust in code
  [ ] Split utils classes and functions  
  [ ] utils.get_bonds use vdw_r in constants  
  [ ] utils.get_angles_and_dihedrals to get Impropers  
  [ ] utils.rot_xyz flag for order of rotations
  [ ] update string function in classes (atom, bond, ...)
  [ ] Check class equality  
  [ ] Molecule ID num in system  
  [ ] system print atoms - to string
  [ ] utils.atom class - print statement - output style
  [ ] Change procrustes name  
  [ ] Merge utils.center_frames into utils.procrustes
  [ ] pretty_xyz -> smooth_frames
  [ ] move utils.colour stuff to another file  
  [ ] sysconst.py should be config file (or merlin to read conf)  
  [ ] separate folder for console commands  
  [ ] separate out optimizers from NEB  
  [ ] install constants and/or paths into merlin  
  [ ] strip out old log file code and implement more generalized one  
  [ ] document LJ, change name (benchmark), move to unit tests  
  [ ] move get lists python files to console commands folder, maybe add to only one file  
  [ ] move gauss_err_fix somewhere (console maybe?) - document
  [ ] overlapping atoms, 180 deg, build linter for g09 and orca
  [ ] OPLS in files to new MD Parameter file - make as class  
  [ ] files.check_net_charge -> utils.py  
  [ ] Lammps.py ?
  [ ] pysub to console commands dir  
  [ ] merge notes into readme documentation  
  [ ] gcube console command  
  [ ] awk console command
3. Divide up documentation between people  
  [ ] James:
  	[ ] Utils classes, utils.rand_rotation, utils.opls stuff
  [ ] Ace:
  	[ ] utils.rot stuff
  [ ] Henry:
  	[ ] utils.improper_angle class distinction
  	[ ] units.py & exceptions
  	[ ] scanDFT, chkDFT  
  	[ ] orca
  	[ ] neb
  	[ ] gauss_err_fix, get_lists
  	[ ] constants  
4. Switch to pep8 style  
  [ ] Fix variables:
  	[ ] Underscore variables and functions
  	[ ] Camal Case for classes
  	[ ] No 'hold' or 'tmp' or 'chk' variables!
  [ ] 160 char line limit  
  [ ] Push variable declaration to top of file/function

5. Subclass Exception to ClancelotException  
6. Sublclass DeprecationWarning to ClancelotDeprecationWarning  
7. Warnings & Exceptions & Deprecation file  

8. How to add new tests (commit hooks)  
9. Pull out NBS  

10. Quantum Espresso Support  
11. LAMMPs Support - abstract out input file  
12. Orca default routes? Best practices  
13. Molecular monte carlo with toehee

14. Glossary  
15. structure and/or parameter (OPLS, ...) database  
