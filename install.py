import os, sys
from getpass import getuser

# In the following list, please ensure you have chosen what you want to install.  By default everything is selected
# If you are installing on a non icse machine, only merlin will install. Recommend turning off all other installations,
# except possibly 'vmd default settings' if you have vmd installed properly

INSTALL_EVERYTHING=False
to_install = {
'vmd':1,
'pysub':1,
'gcube':1,
'jsub':1,
'jdel':1,
'chkDFT':1,
'scanDFT':1,
'junest (formerly juju)':0,
'anaconda':0, 				# a Python 2.7.9 distribution that installs to ~/anaconda
'vmd default settings':0,	# improves the default settings of vmd
'file_browser':1, 			# set the file browser not to open a new window per folder
'merlin':1,
'sublime_text_3_build_3083':0,
'prnt':1 # ONLY install if lpstat -p -d returns no available printers
}

####################################################################################################################
if INSTALL_EVERYTHING:
	for x in to_install:
		to_install[x] = 1

if to_install['file_browser']:
	os.system('gconftool-2   --type bool --set /apps/nautilus/preferences/always_use_browser true')

# Use HOMEDIR instead of /fs/home/USERNAME in order to install on non icse machines
USERNAME = getuser()
HOMEDIR = os.path.expanduser('~')

INSTALLDIR = os.getcwd()
if INSTALLDIR[-1] != '/': INSTALLDIR += '/' # Ensure there is a trailing slash
ZSHRC = HOMEDIR+'/.zshrc'
ZSH_CLANCELOT = HOMEDIR + '/.zsh_clancelot'
BASHRC = HOMEDIR+'/.bashrc'
temp=open(ZSHRC)
zshrc_string=temp.read()
temp.close()

def zshrc_check_add(st,path,zshrc_contents):
	if st not in zshrc_contents:
		f=open(ZSHRC,'a')
		f.write('\n'+st+'\n')
		f.close

if 'source ~/.zsh_clancelot' not in zshrc_string:
	f=open(ZSHRC,'a')
	f.write('''\n# The following loads the Clancelot Config File
if [ -f ~/.zsh_clancelot ]; then
    source ~/.zsh_clancelot
else
    print '404: ~/.zsh_clancelot not found.'
fi
''')
	f.close()
	f=open(BASHRC,'a')
	f.write('''\n# The following loads the Clancelot Config File
if [ -f ~/.zsh_clancelot ]; then
    source ~/.zsh_clancelot
else
    print '404: ~/.zsh_clancelot not found.'
fi
''')
	f.close()


os.system('mkdir -p '+INSTALLDIR) # Ensure the install directory is made
for key in to_install: # Make directories for what we want to install
	if key == 'chkDFT': continue
	if key == 'jdel': continue
	if key == 'jsub': continue
	if key == 'vmd': continue
	if key == 'scanDFT': continue
	if key == 'junest (formerly juju)': continue
	if key == 'python 2.7.10': continue
	if key == 'numpy': continue
	if key == 'scipy': continue
	if key == 'cython-0.22': continue
	if key == 'vmd default settings': continue
	if key == 'file_browser': continue
	if key == 'matplotlib': continue
	if key == 'anaconda': continue
	if key == 'sublime_text_3_build_3083': continue
	if key == 'merlin': continue
	if key == 'prnt': continue
	if to_install[key]: os.system('mkdir -p '+INSTALLDIR+key+'/')


# Give Bash Capabilities to ZSH
f = open(ZSH_CLANCELOT,'w+')
f.write('''###############################################################
############### THE FOLLOWING IS FOR CLANCELOT ################
###############################################################
# Append Path
export PYTHONPATH=$$$$$$/tools:$PYTHONPATH
export PATH=/fs/europa/g_pc/orca_3_0_3_linux_x86-64/:$PATH

# Aliases for our tools
alias get_ext_list='python $$$$$$/tools/get_ext_list.py'
alias get_gauss_list='python $$$$$$/tools/get_gauss_list.py'
alias get_orca_list='python $$$$$$/tools/get_orca_list.py'
alias get_jlist='python $$$$$$/tools/get_jlist.py'

# Bind keys
bindkey '^[[3~' delete-char
bindkey '^[OH' beginning-of-line
bindkey '^[OF' end-of-line
bindkey ';5C' emacs-forward-word
bindkey ';5D' emacs-backward-word

'''.replace('$$$$$$/',INSTALLDIR))

f.close()
f = open(ZSH_CLANCELOT,'a')
f.write('''autoload bashcompinit # Let us use bash commands
bashcompinit

# Note, we will provide an explanation of how these functions work using _jAutoTab, but the rest
# of these functions will be largely uncommented

# This function provides auto-tab for the jlist
_jAutoTab() # By convention, the function name starts with an underscore.
{
local cur # Pointer to current completion word.
local JLIST # Pointer to a variable that will hold your list
COMPREPLY=() # Array variable storing the possible completions.
cur=${COMP_WORDS[COMP_CWORD]}
# Note, you can get the list however you choose. In this example, we call an aliased python file and pass it
# the current working directory so that it can return a list. These lists are in the format of space
# separated strings. For example: 'file1 file2 file3'. Note, we just need to print the list to screen from
# this python file, not return it.
JLIST=$(get_jlist)
# This is the function that will determine what words from your JLIST will be displayed for autocomplete
case "$cur" in
*)
COMPREPLY=( $( compgen -W '$JLIST' $cur ) );; # You need to enter your list here
esac
return 0
}

_pyAutoTab()
{
local cur 
local PYLIST 
COMPREPLY=() 
cur=${COMP_WORDS[COMP_CWORD]}
PYLIST=$(get_ext_list .py $PWD)
case "$cur" in
*)
COMPREPLY=( $( compgen -W '$PYLIST' $cur ) );;
esac
return 0
}

_nbsAutoTab()
{
local cur 
local NBSLIST 
COMPREPLY=() 
cur=${COMP_WORDS[COMP_CWORD]}
NBSLIST=$(get_ext_list .nbs $PWD)
case "$cur" in
*)
COMPREPLY=( $( compgen -W '$NBSLIST' $cur ) );;
esac
return 0
}

_gaussAutoTab()
{
local cur 
local GAUSSLIST 
COMPREPLY=() 
cur=${COMP_WORDS[COMP_CWORD]}
GAUSSLIST=$(get_gauss_list $PWD)
case "$cur" in
*)
COMPREPLY=( $( compgen -W '$GAUSSLIST' $cur ) );;
esac
return 0
}

_orcaAutoTab()
{
local cur 
local ORCALIST 
COMPREPLY=() 
cur=${COMP_WORDS[COMP_CWORD]}
ORCALIST=$(get_orca_list $PWD)
case "$cur" in
*)
COMPREPLY=( $( compgen -W '$ORCALIST' $cur ) );;
esac
return 0
}''')
f.close


f = open(ZSH_CLANCELOT,'a')
f.write('\n\n')

if to_install['vmd']: f.write("alias vmd='/fs/europa/g_pc/vmd/bin/vmd'\n\n")
if to_install['pysub']:
	f.write("alias pysub='"+INSTALLDIR+"pysub/pysub.sh'\n")
	f.write('complete -F _pyAutoTab '+INSTALLDIR+'pysub/pysub.sh\n\n')
	g = open(INSTALLDIR+'pysub/pysub.sh','w')
	g.write('python '+INSTALLDIR+'''pysub/pysub.py $PWD'/' $@''')
	g.close()
	g = open(INSTALLDIR+'pysub/pysub.py','w')
	g.write('''import os, sys
from time import sleep

s_hold = \'\'\'#!/bin/bash
##NBS-name: "$$$$$$"
##NBS-nproc: 1
##NBS-queue: "batch"

rm ^^^^^^^$$$$$$.log

source '''+HOMEDIR+'''/.zshrc

'''+HOMEDIR+'''/anaconda/bin/python2.7 -u ^^^^^^^$$$$$$.py >> ^^^^^^^$$$$$$.log 2>&1
\'\'\'
for i,s in enumerate(sys.argv):
	if i == 0: continue
	if i == 1: continue
	f=open('pysub.nbs','w')
	f.write(s_hold.replace('$$$$$$',s[:-3]).replace('^^^^^^^',sys.argv[1][:sys.argv[1].rfind('/')+1]))
	f.close()
	os.system('jsub pysub.nbs')

os.system('rm pysub.nbs')''')
	g.close()
	os.system('chmod 755 '+INSTALLDIR+'pysub/pysub.sh')
if to_install['gcube']:
	f.write("alias gcube='"+INSTALLDIR+"gcube/cube.sh'\n")
	f.write('complete -F _gaussAutoTab '+INSTALLDIR+'gcube/cube.sh\n\n')
	g = open(INSTALLDIR+'gcube/cube.sh','w')
	g.write('python '+INSTALLDIR+'''gcube/cube.py $PWD'/' $@''')
	g.close()
	g = open(INSTALLDIR+'gcube/cube.py','w')
	g.write('''from merlin import *
from subprocess import Popen

old_job = sys.argv[2]

if not os.path.exists('gaussian/%s.chk' % old_job):
	print 'Fatal error: file "gaussian/%s.chk" does not exist.' % old_job
	exit()

if not g09.parse_atoms(old_job):
	print 'Fatal error: "%s" is not converged. gcube does not work on unconverged jobs.' % old_job
	exit()

# Get the file to check
if not os.path.exists('gaussian/%s.fchk' % old_job):
	print 'Making gaussian/%s.fchk' % old_job
	Popen('/usr/local/gaussian/g09/g09/formchk gaussian/%s.chk gaussian/%s.fchk' % (old_job,old_job), shell=True).wait()

# Make the density and potential cube files
if not os.path.exists('gaussian/%s.cube' % (old_job+'_d')):
	print 'Making gaussian/%s.cube' % (old_job+'_d')
	Popen('/usr/local/gaussian/g09/g09/cubegen 0 density gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_d'), shell=True).wait()

if not os.path.exists('gaussian/%s.cube' % (old_job+'_p')):
	print 'Making gaussian/%s.cube' % (old_job+'_p')
	Popen('/usr/local/gaussian/g09/g09/cubegen 0 potential gaussian/%s.fchk gaussian/%s.cube 0 h'% (old_job,old_job+'_p'), shell=True).wait()

if not (os.path.exists('gaussian/%s.cube' % (old_job+'_p')) and os.path.exists('gaussian/%s.cube' % (old_job+'_d')) ):
	print 'Fatal error: cube files not created'
	exit()

vmd_file = \'\'\'# Type logfile console into console to see all commands

# Get data
mol new gaussian/$$FPTR$$_d.cube
mol addfile gaussian/$$FPTR$$_p.cube

# Adjust first rep
mol modcolor 0 0 element
mol modstyle 0 0 CPK

# Adjust second rep
mol addrep 0
mol modcolor 1 0 Volume 1
mol modstyle 1 0 Isosurface 0.040000 0 0 0 1 1
mol modmaterial 1 0 Transparent\'\'\'.replace('$$FPTR$$',old_job)

f = open('tmp.vmd','w')
f.write(vmd_file)
f.close()

Popen('/fs/europa/g_pc/vmd-1.9 -e tmp.vmd', shell=True)\n''')
	g.close()
	os.system('chmod 755 '+INSTALLDIR+'gcube/cube.sh')
if to_install['jsub']: f.write('complete -F _nbsAutoTab jsub\n\n')
if to_install['jdel']: f.write('complete -F _jAutoTab jdel\n\n')
if to_install['chkDFT']:
	f.write("alias chkDFT='python "+INSTALLDIR+"console_scripts/chkDFT.py'\n")
	f.write('''alias viewg='function _viewg(){chkDFT $1 -dft g09 -v $@};_viewg'
alias viewo='function _viewo(){chkDFT $1 -dft orca -v $@};_viewo'
alias chkg='function _chkg(){chkDFT $1 -dft g09 $@};_chkg'
alias ggedit='function _ggedit(){gedit orca/$1/$1.log &};_ggedit'
alias geditg='ggedit'
alias gtail='function _gtail(){tail orca/$1/$1.log $2 $3};_gtail'
alias tailg='gtail'
alias chko='function _chko(){chkDFT $1 -dft orca $@};_chko'
alias ogedit='function _ogedit(){gedit orca/$1/$1.out &};_ogedit'
alias gedito='ogedit'
alias otail='function _otail(){tail orca/$1/$1.out $2 $3};_otail'
alias tailo='otail'\n''')
if to_install['merlin']: f.write("alias merlin='python -i "+INSTALLDIR+"tools/merlin.py'\n")
if to_install['scanDFT']: f.write("\nalias scanDFT='python "+INSTALLDIR+"console_scripts/scanDFT.py'\n")
if to_install['prnt']: f.write('''alias prnt='function _prnt(){ssh asimov "lpr -P hplj4525-365 -o sides=two-sided-long-edge -o InputSlot=Tray2 $PWD/$1;logout";echo "Printed..."};_prnt'\n''')
f.write('''\n###############################################################
################## END OF THE CLANCELOT CODE ##################
###############################################################''')
f.close()

downloaded_tarball=False

#####The following is modified from http://code.activestate.com/recipes/134892/
#####As opposed to raw_input(), getchar() does not wait for a new line (*enter*)
def getchar():
    import tty, termios
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        tty.setraw(sys.stdin.fileno())
        ch = sys.stdin.read(1)
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
    return ch
#####


def reinstall(str):
	while True:
		print(str+' y/n:')
		resp=getchar()
		if resp=='y':
			return True
		if resp=='n':
			return False

if to_install['vmd default settings']:
	if os.path.exists(HOMEDIR+'/.vmdrc'):
		if reinstall('Previous vmd settings found, backup old ~/.vmdrc to ~/.vmdrc_history?'):
			os.system('mv ~/.vmdrc ~/.vmdrc_history')
	os.system('cp vmdrc_default.txt ~/.vmdrc')

def anaconda_install():
	os.system('wget -P ~/lib/ https://repo.continuum.io/archive/Anaconda-2.2.0-Linux-x86_64.sh')
	os.system('bash ~/lib/Anaconda-2.2.0-Linux-x86_64.sh -fb')
	zshrc_check_add('export PATH=~/anaconda/bin:$PATH',ZSHRC,zshrc_string)
	#zshrc_check_add("export PYTHONPATH='"+INSTALLDIR+"/tools/'",ZSHRC,zshrc_string)
	os.system('rm ~/lib/Anaconda-2.2.0-Linux-x86_64.sh')

def sublime_install():
	os.system('mkdir -p '+HOMEDIR+'/lib')
	os.system('wget -P ~/lib/ http://c758482.r82.cf2.rackcdn.com/sublime_text_3_build_3083_x64.tar.bz2')
	os.system('tar xvf '+HOMEDIR+'/lib/sublime_text_3_build_3083_x64.tar.bz2 -C '+HOMEDIR+'/lib/')
	zshrc_check_add("alias sublime='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)
	zshrc_check_add("alias subl='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)

def junest_install():
	os.system('git clone git://github.com/fsquillace/juju ~/juju --quiet')
	zshrc_check_add("export PATH=~/juju/bin:$PATH",ZSHRC,zshrc_string)
	zshrc_check_add("alias juju='junest'",ZSHRC,zshrc_string)

if to_install['anaconda']:
	if os.path.exists(HOMEDIR+'/anaconda') and os.path.isdir(HOMEDIR+'/anaconda'):
		if reinstall('Previous installation found, reinstall Anaconda (Python 2.7.9 and packages)?'):
			anaconda_install()
		else:
			print('...SKIPPING ANACONDA (RE)INSTALLATION...')
	else:
		anaconda_install()
else:
	if os.path.exists(HOMEDIR+'/anaconda') and os.path.isdir(HOMEDIR+'/anaconda'):
		zshrc_check_add('export PATH=~/anaconda/bin:$PATH',ZSHRC,zshrc_string)


if to_install['sublime_text_3_build_3083']:
	if os.path.exists(HOMEDIR+'/lib/sublime_text_3') and os.path.isdir(HOMEDIR+'/lib/sublime_text_3'):
		if reinstall('Previous installation found, reinstall Sublime Text 3?'):
			os.system('rm -rf '+HOMEDIR+'/lib/sublime_text_3')
			sublime_install()
			downloaded_tarball=True
		else:
			print('...SKIPPING SUBLIME (RE)INSTALLATION...')
	else:
		sublime_install()
		downloaded_tarball=True
else:
	if os.path.exists(HOMEDIR+'/lib/sublime_text_3') and os.path.isdir(HOMEDIR+'/lib/sublime_text_3'):
		zshrc_check_add("alias sublime='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)
		zshrc_check_add("alias subl='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)

if to_install['junest (formerly juju)']: 
	if os.path.exists(HOMEDIR+'/juju') and os.path.isdir(HOMEDIR+'/juju'):
		if reinstall('Previous installation found, reinstall juju/junest?'):
			os.system('mv '+HOMEDIR+'/juju '+HOMEDIR+'/.trash/')
			if os.path.exists(HOMEDIR+'/.junest') and os.path.isdir(HOMEDIR+'/.junest'):
				os.system('mv '+HOMEDIR+'/.junest '+HOMEDIR+'/.trash/')
			junest_install()
			print("\nTo finish installing 'junest' please run:\n'pacman -Syyu pacman-mirrorlist && pacman -S gtk2 avogadro grep make ttf-liberation gedit'\n\n(when prompted for GL version, pick option 2, nvidia)\n\n\n")
			os.system("zsh -c 'junest -f'")
		else:
			print('...SKIPPING JUJU/JUNEST (RE)INSTALLATION...')
	else:
		junest_install()
		print("\nTo finish installing 'junest' please run:\n'pacman -Syyu pacman-mirrorlist && pacman -S gtk2 avogadro grep make ttf-liberation gedit'\n\n(when prompted for GL version, pick option 2, nvidia)\n\n\n")
		os.system("zsh -c 'junest -f'")
else:
	if os.path.exists(HOMEDIR+'/juju') and os.path.isdir(HOMEDIR+'/juju'):
		zshrc_check_add("export PATH=~/juju/bin:$PATH",ZSHRC,zshrc_string)
		zshrc_check_add("alias juju='junest'",ZSHRC,zshrc_string)

if downloaded_tarball:
	print('Removing previously downloaded tarballs')
	os.system('rm -i '+HOMEDIR+'/lib/*.tar.*')

os.system('cp pre-commit.sh .git/hooks/pre-commit') #copy pre-commit hook script to expected location
os.system('chmod +x .git/hooks/pre-commit') #make pre-commit hook script executable

print("\n\n--------------Installation Finished--------------\nPlease reopen Terminal to apply changes.")

