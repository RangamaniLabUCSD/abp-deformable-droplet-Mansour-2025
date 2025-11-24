import Traj_analyze_fracoccupied as main
import sys
N = int(sys.argv[1])
Rval = int(sys.argv[2])
Nreps = int(sys.argv[3])
print(len(sys.argv),flush=True)
if(len(sys.argv)==5):
    main.frontend_filprops_Rval(N, Rval, Nreps, sys.argv[4])
elif(len(sys.argv)==6):
    main.frontend_filprops_Rval(N, Rval, Nreps, sys.argv[4], sys.argv[5])
else:
    main.frontend_filprops_Rval(N, Rval, Nreps)