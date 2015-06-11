import os, sys
import math

# These imports are for the notify and watch functions -----
from subprocess import Popen, PIPE
from time import sleep
import smtplib
from email.mime.text import MIMEText
# ----------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# These two functions help notify me when a simulation job has finished
def notify(job_name, cell):
	msg = MIMEText(job_name + " has finished")

	msg['Subject'] = job_name + " has finished"
	msg['From'] = "g_pc@cbe.cornell.edu"
	msg['To'] = cell

	# Send the message via our own SMTP server, but don't include the envelope header.
	s = smtplib.SMTP('localhost')
	s.sendmail("g_pc@cbe.cornell.edu", [cell], msg.as_string())
	s.quit()

def watch(job_name, cell):
	while 1:
		p = Popen(['jlist'], stdout=PIPE)
		s = p.stdout.read().split()
		dne = True
		for i in range(len(s)):
			if((s[i] == job_name) and ((i + 2) < len(s)) and (s[i+2] == 'RUNNING')):
				dne = False
		if(dne):
			break
		sleep(10)

	notify(job_name, cell)
# --------------------------------------------------------------------------------------------------------------
# Simple checks on numeric status of string
def isInt(s):
    try: 
        int(s)
        return True
    except ValueError: return False

def isNum(s):
    try: 
        float(s)
        return True
    except ValueError: return False
# --------------------------------------------------------------------------------------------------------------