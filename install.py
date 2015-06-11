import os, sys
from getpass import getuser

# In the following list, please ensure you have chosen what you want to install.  By default everything is selected
to_install = {
'vmd':1,
'pysub':1,
'jsub':1,
'jdel':1,
'viewg':1,
'chkg':1,
'chkg_all':1,
'scang':1,
'junest (formerly juju)':1,
'anaconda':1, 				# a Python 2.7.9 distribution that installs to ~/anaconda
'vmd default settings':1,	# improves the default settings of vmd
'file_browser':1, 			# set the file browser not to open a new window per folder
'merlin':1,
'sublime_text_3_build_3083':1
}

####################################################################################################################
if to_install['file_browser']:
	os.system('gconftool-2   --type bool --set /apps/nautilus/preferences/always_use_browser true')

USERNAME = getuser()


INSTALLDIR = os.getcwd()
if INSTALLDIR[-1] != '/': INSTALLDIR += '/' # Ensure there is a trailing slash
ZSHRC = '/fs/home/'+USERNAME+'/.zshrc'
ZSH_CLANCELOT = '/fs/home/' + USERNAME + '/.zsh_clancelot'
BASHRC = '/fs/home/'+USERNAME+'/.bashrc'
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
	if key == 'chkg_all': continue
	if key == 'jdel': continue
	if key == 'jsub': continue
	if key == 'vmd': continue
	if key == 'scang': continue
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
	if to_install[key]: os.system('mkdir -p '+INSTALLDIR+key+'/')


# Give Bash Capabilities to ZSH
f = open(ZSH_CLANCELOT,'w+')
f.write('''###############################################################
############### THE FOLLOWING IS FOR CLANCELOT ################
###############################################################
# Append Path
export PYTHONPATH=$$$$$$/tools:$PYTHONPATH

# Aliases for our tools
alias get_ext_list='python $$$$$$/tools/get_ext_list.py'
alias get_gauss_list='python $$$$$$/tools/get_gauss_list.py'
alias get_jlist='python $$$$$$/tools/get_jlist.py'
alias merlin='python -i $$$$$$/tools/merlin.py'

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
}''')
f.close


f = open(ZSH_CLANCELOT,'a')
f.write('\n\n')

if to_install['vmd']: f.write("alias vmd='/fs/europa/g_pc/vmd-1.9'\n\n")
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

python ^^^^^^^$$$$$$.py
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
if to_install['jsub']: f.write('complete -F _nbsAutoTab jsub\n\n')
if to_install['jdel']: f.write('complete -F _jAutoTab jdel\n\n')
if to_install['viewg']:
	f.write("alias viewg='"+INSTALLDIR+"viewg/viewg.sh'\n")
	f.write('complete -F _gaussAutoTab '+INSTALLDIR+'viewg/viewg.sh\n\n')
	g = open(INSTALLDIR+'viewg/viewg.sh','w')
	g.write('python '+INSTALLDIR+'tools/view_gaussian.py $@')
	g.close()
	os.system('chmod 755 '+INSTALLDIR+'viewg/viewg.sh')
if to_install['chkg']:
	f.write("alias chkg='"+INSTALLDIR+"chkg/chkg.sh'\n")
	f.write('complete -F _gaussAutoTab '+INSTALLDIR+'chkg/chkg.sh\n\n')
	g = open(INSTALLDIR+'chkg/chkg.sh','w')
	g.write('python '+INSTALLDIR+'tools/chkg.py $@')
	g.close()
	os.system('chmod 755 '+INSTALLDIR+'chkg/chkg.sh')
if to_install['chkg_all']: f.write("alias chkg_all='python "+INSTALLDIR+"tools/chkg_all.py'\n")
if to_install['merlin']: f.write("alias merlin='python "+INSTALLDIR+"tools/merlin.py'\n")
if to_install['scang']: f.write("\nalias scang='python /fs/home/jms875/scan.py'\n")

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
	if os.path.exists('/fs/home/'+USERNAME+'/.vmdrc'):
		if reinstall('Previous vmd settings found, backup old ~/.vmdrc to ~/.vmdrc_history?'):
			os.system('mv ~/.vmdrc ~/.vmdrc_history')
	os.system('cp vmdrc_default.txt ~/.vmdrc')

def anaconda_install():
	os.system('wget -P ~/lib/ https://repo.continuum.io/archive/Anaconda-2.2.0-Linux-x86_64.sh')
	os.system('bash ~/lib/Anaconda-2.2.0-Linux-x86_64.sh -fb')
	zshrc_check_add('export PATH=~/anaconda/bin:$PATH',ZSHRC,zshrc_string)
	zshrc_check_add("export PYTHONPATH=''",ZSHRC,zshrc_string)
	os.system('rm ~/lib/Anaconda-2.2.0-Linux-x86_64.sh')

def sublime_install():
	os.system('mkdir -p /fs/home/'+USERNAME+'/lib')
	os.system('wget -P ~/lib/ http://c758482.r82.cf2.rackcdn.com/sublime_text_3_build_3083_x64.tar.bz2')
	os.system('tar xvf /fs/home/'+USERNAME+'/lib/sublime_text_3_build_3083_x64.tar.bz2 -C /fs/home/' + USERNAME + '/lib/')
	zshrc_check_add("alias sublime='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)
	zshrc_check_add("alias subl='~/lib/sublime_text_3/sublime_text",ZSHRC,zshrc_string)

def junest_install():
	os.system('git clone git://github.com/fsquillace/juju ~/juju --quiet')
	zshrc_check_add("export PATH=~/juju/bin:$PATH",ZSHRC,zshrc_string)
	zshrc_check_add("alias juju='junest'",ZSHRC,zshrc_string)

if to_install['anaconda']:
	if os.path.exists('/fs/home/'+USERNAME+'/anaconda') & os.path.isdir('/fs/home/'+USERNAME+'/anaconda'):
		if reinstall('Previous installation found, reinstall Anaconda (Python 2.7.9 and packages)?'):
			anaconda_install()
		else:
			print('...SKIPPING ANACONDA (RE)INSTALLATION...')
	else:
		anaconda_install()

if to_install['sublime_text_3_build_3083']:
	if os.path.exists('/fs/home/'+USERNAME+'/lib/sublime_text_3') & os.path.isdir('/fs/home/'+USERNAME+'/lib/sublime_text_3'):
		if reinstall('Previous installation found, reinstall Sublime Text 3?'):
			os.system('rm -rf /fs/home/'+USERNAME+'/lib/sublime_text_3')
			sublime_install()
			downloaded_tarball=True
		else:
			print('...SKIPPING SUBLIME (RE)INSTALLATION...')
	else:
		sublime_install()
		downloaded_tarball=True

if to_install['junest (formerly juju)']: 
	if os.path.exists('/fs/home/'+USERNAME+'/juju') & os.path.isdir('/fs/home/'+USERNAME+'/juju'):
		if reinstall('Previous installation found, reinstall juju/junest?'):
			os.system('mv /fs/home/'+USERNAME+'/juju /fs/home/'+USERNAME+'/.trash/')
			if os.path.exists('/fs/home/'+USERNAME+'/.junest') & os.path.isdir('/fs/home/'+USERNAME+'/.junest'):
				os.system('mv /fs/home/'+USERNAME+'/.junest /fs/home/'+USERNAME+'/.trash/')
			junest_install()
			print("\nTo finish installing 'junest' please run:\n'pacman -Syyu pacman-mirrorlist && pacman -S gtk2 avogadro grep make ttf-liberation gedit'\n\n(when prompted for GL version, pick option 2, nvidia)\n\n\n")
			os.system("zsh -c 'junest -f'")
		else:
			print('...SKIPPING JUJU/JUNEST (RE)INSTALLATION...')
	else:
		junest_install()
		print("\nTo finish installing 'junest' please run:\n'pacman -Syyu pacman-mirrorlist && pacman -S gtk2 avogadro grep make ttf-liberation gedit'\n\n(when prompted for GL version, pick option 2, nvidia)\n\n\n")
		os.system("zsh -c 'junest -f'")

if downloaded_tarball:
	print('Removing previously downloaded tarballs')
	os.system('rm -i /fs/home/'+USERNAME+'/lib/*.tar.*')

print("\n\n--------------Installation Finished--------------\nPlease reopen Terminal to apply changes.")

