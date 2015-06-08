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
'python 2.7.10':1,
'cython-0.22':1,
'numpy':1,
'scipy':1,
'vmd default settings':1
}
# Is this your first time running this script? (To avoid redundant additions to .zshrc)
first_time = 1

####################################################################################################################
USERNAME = getuser()


INSTALLDIR = os.getcwd()
if INSTALLDIR[-1] != '/': INSTALLDIR += '/' # Ensure there is a trailing slash
ZSHRC = '/fs/home/'+USERNAME+'/.zshrc'
ZSH_CLANCELOT = '/fs/home/' + USERNAME + '/.zsh_clancelot'
BASHRC = '/fs/home/'+USERNAME+'/.bashrc'

if first_time:
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
	if to_install[key]: os.system('mkdir -p '+INSTALLDIR+key+'/')


# Give Bash Capabilities to ZSH
f = open(ZSH_CLANCELOT,'w+')
f.write('''###############################################################
############### THE FOLLOWING IS FOR CLANCELOT ################
###############################################################

# Aliases for our tools
alias get_ext_list='python $$$$$$/tools/get_ext_list.py'
alias get_gauss_list='python $$$$$$/tools/get_gauss_list.py'
alias get_jlist='python $$$$$$/tools/get_jlist.py'

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
if to_install['scang']: f.write("\nalias scang='python /fs/home/jms875/scan.py'\n")

f.write('''\n###############################################################
################## END OF THE CLANCELOT CODE ##################
###############################################################''')
f.close()

if to_install['vmd default settings']:
	os.system('cp vmdrc_default.txt ~/.vmdrc')

if to_install['python 2.7.10']:
	os.system('mkdir -p /fs/home/' + USERNAME + '/lib')
	os.system('wget -P /fs/home/'+USERNAME+'/lib/ https://www.python.org/ftp/python/2.7.10/Python-2.7.10.tar.xz')
	os.system('tar xvf /fs/home/'+USERNAME+'/lib/Python-2.7.10.tar.xz -C /fs/home/' + USERNAME + '/lib/')
	os.system('cd ~/lib/Python-2.7.10/')
	os.chdir('/fs/home/'+USERNAME+'/lib/Python-2.7.10/')
	os.system('bash /fs/home/'+USERNAME+'/lib/Python-2.7.10/configure --prefix=/fs/home/'+USERNAME+'/lib/Python-2.7.10/')
	os.system('make')
	os.system('make install')
	f = open(ZSHRC,'a')
	f.write("\nexport PATH=$HOME/lib/Python-2.7.10/bin:$PATH\n")
	f.close()
	os.system('export PATH=~/lib/Python-2.7.10/bin:$PATH')
	
if to_install['cython-0.22']:
	os.system('mkdir -p /fs/home/' + USERNAME + '/lib')
	os.system('wget -P ~/lib/ http://cython.org/release/Cython-0.22.tar.gz')
	os.system('tar xvf /fs/home/'+USERNAME+'/lib/Cython-0.22.tar.gz -C /fs/home/' + USERNAME + '/lib/')
	os.system('cd ~/lib/Cython-0.22')
	os.chdir('/fs/home/'+USERNAME+'/lib/Cython-0.22')
	os.system('python setup.py install')
	f = open(ZSHRC,'a')
	f.write("\nexport PATH=~/lib/Cython-0.22/bin:$PATH\n")
	f.close()
	os.system('export PATH=~/lib/Cython-0.22/bin:$PATH')

if to_install['numpy']:
	os.system('mkdir -p /fs/home/' + USERNAME + '/lib/Python-2.7.10')
	os.system('git clone git://github.com/numpy/numpy.git ~/lib/numpy')
	os.chdir('/fs/home/'+USERNAME+'/lib/numpy')
	os.system('python setup.py install --home=~/lib/Python-2.7.10/')
	os.system('export PYTHONPATH=~/lib/Python-2.7.10/lib/python:$PYTHONPATH')
	f = open(ZSHRC,'a')
	f.write("\nexport PYTHONPATH=~/lib/Python-2.7.10/lib/python:$PYTHONPATH\n")
	f.close()
	
if to_install['scipy']:
	os.system('mkdir -p /fs/home/' + USERNAME + '/lib/Python-2.7.10')
	os.system('git clone git://github.com/scipy/scipy.git ~/lib/scipy')
	os.chdir('/fs/home/'+USERNAME+'/lib/scipy')
	os.system('python setup.py install --home=~/lib/Python-2.7.10/')
	os.system('export PYTHONPATH=~/lib/Python-2.7.10/lib/python:$PYTHONPATH')
	f = open(ZSHRC,'a')
	f.write("\nexport PYTHONPATH=~/lib/Python-2.7.10/lib/python:$PYTHONPATH\n")
	f.close()

if to_install['junest (formerly juju)']: 
	os.system('git clone git://github.com/fsquillace/juju ~/juju --quiet')
	f = open(ZSHRC,'a')
	f.write("\nexport PATH=~/juju/bin:$PATH\n")
#	f.write("alias juju='junest'")
	f.close()
	print("\nTo finish installing 'junest' please run:\n'pacman -Syyu pacman-mirrorlist && pacman -S gtk2 avogadro grep make ttf-liberation gedit'\n\n(when prompted for GL version, pick option 2, nvidia)\n\n\n")
	os.system("zsh -c 'junest -f'")


