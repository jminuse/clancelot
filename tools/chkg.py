import sys, log

if (len(sys.argv)==2): log.chkg_E(sys.argv[1],record=False)
if (len(sys.argv)==3): log.chkg_E(sys.argv[1],unit=sys.argv[2],record=False)
if (len(sys.argv)==4): log.chkg_E(sys.argv[1],unit=sys.argv[2],record=bool(int(sys.argv[3])))
if (len(sys.argv)>4): print('Too many arguments')