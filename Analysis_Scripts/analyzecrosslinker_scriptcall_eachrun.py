import Traj_analyze_fracoccupied as main
import sys
N = int(sys.argv[1])
Rval = int(sys.argv[2])
repid = int(sys.argv[3])
print(len(sys.argv),flush=True)
if(len(sys.argv)==5):
    main.frontend_cross_set_Rval_repid(N, Rval, repid, sys.argv[4])
elif(len(sys.argv)==6):
    main.frontend_cross_set_Rval_repid(N, Rval, repid, sys.argv[4], sys.argv[5])
else:
    main.frontend_cross_set_Rval_repid(N, Rval, repid)
