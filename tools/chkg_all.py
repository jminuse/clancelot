import sys
from log import chkg_E_all

if (len(sys.argv)==1): chkg_E_all()
if (len(sys.argv)==2): chkg_E_all(u=sys.argv[1])
