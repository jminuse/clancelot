from merlin import *
import re


def spaced_print(sOut,delim=['\t',' '],buf=4):
        s_len = []
        if type(sOut) == str: sOut = sOut.split('\n')
        if type(delim) == list: delim = ''.join([d+'|' for d in delim])[:-1]
        # Get the longest length in the column
        for i,s in enumerate(sOut):
                s = re.split(delim,s)
                for j,ss in enumerate(s):
                        try: s_len[j] = len(ss) if len(ss)>s_len[j] else s_len[j] # This makes the part of the list for each column the longest length
                        except: s_len.append(len(ss)) # If we are creating a new column this happens
        for i in range(len(s_len)): s_len[i] += buf # Now we add a buffer to each column

        # Compile string output
        for i,s in enumerate(sOut):
                s = re.split(delim,s)
                for j,ss in enumerate(s): s[j] = ss + ''.join([' ']*(s_len[j]-len(ss)))
                sOut[i] = ''.join(s)

        return '\n'.join(sOut)

if len(sys.argv) < 2:
	print("Not enough arguments. Input either 1 indicating the simulation or 2 indicating simulation and units.")
	sys.exit()

a,e,t,conv = orca.parse_atoms(sys.argv[1],get_atoms=True,get_energy=True,get_charges=False,get_time=True,check_convergence=True,parse_all=True)

u = 'Ha'
if len(sys.argv) > 2: u2 = sys.argv[2]
else: u2 = 'Ha'

head = 'Job Name: %s\n' % sys.argv[1]
head += 'Energy Data Points: %d\n' % len(e)
Ener = str(units.convert_energy(u, u2, e[-2] - e[-3]))
head += 'dE 2nd last = %s %s\n' % (Ener,u2)
Ener = str(units.convert_energy(u, u2, e[-1] - e[-2]))
head += 'dE last = %s %s\n' % (Ener,u2)
Ener = str(units.convert_energy(u, u2, e[-1]))
head += 'Last Energy = %s %s' % (Ener,u2)
body = ''
for c in conv[-1]: body += '%s\t%s\t%s\t%s\n' % (c[0],str(c[1]),str(c[2]),str(c[3]))
body = body[:-1]
body = spaced_print(body, delim='\t')

s = '\n'+''.join(['-']*len(body.split('\n')[0]))+'\n'

print(s+head+s+body+s)
