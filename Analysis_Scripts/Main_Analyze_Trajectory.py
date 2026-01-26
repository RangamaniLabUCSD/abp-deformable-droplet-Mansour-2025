######## This python scipt helps analyze the trajectories from Cytosim

######## Imports ########

import math
import numpy as np
import os
import pandas as pd
import icosphere as ic
from datetime import datetime
from collections import Counter
from collections import defaultdict
from joblib import Parallel, delayed
from scipy.stats import gaussian_kde
import subprocess

######## Classes ########

class replicate:
    ridx = 0;
    Nsnaps = 0;
    snap = [];
    timevector = []
    def __init__(self, rid):
        self.ridx = rid
        self.Nsnaps = 0
        self.snap = [];
        self.timevector =[];
    def addsnapshot(self):
        self.snap.append(snapshot(1))

class snapshot:
    sidx = 0;
    spacecoord=[]
    filcoord=[]
    crossboundcoord=[];
    crossboundblobid = [];
    handid = [];
    dimerid = [];
    boundid = [];
    partnerid = [];
    crossboundfilid = [];
    solidcoord = [];
    selfblobid = []
    def __init__(self, sid):
        self.sidx = sid
        self.spacecoord = [];
        self.filcoord = [];
        self.crossboundcoord = [];
        self.crossboundblobid = [];
        self.handid = [];
        self.dimerid = [];
        self.boundid = [];
        self.partnerid = [];
        self.crossboundfilid = [];
        self.solidcoord = [];
        self.selfblobid = [];

######## Functions ########

def readcytosimtraj (*args):
    tag = ''
    filepath = args[0]
    filename = 'objects.cmo'
    if(len(args)==2):
        tag= args[1]
    # Get number of solids
    test = subprocess.run(["bash", "./get_nsolids.sh", "-f", filepath], capture_output=True)
    nsolids = int(test.stdout)
    print("Nsolids in trajectory = "+str(nsolids), flush=True)
    r=[];
    Nruns =1 
    for  i in range(0,Nruns):
        r.append(replicate(i))
    printstatus = False; 
    recordStatus = False; 
    ridx = 0;
    sidx=-1;
    # Open file
    fptr = open(filepath+filename,'r', encoding="utf8", errors="ignore")
    for line in fptr:
        if(printstatus):
            print(line)
        # Simulation snapshot time
        if('time' in line):
            t = float((line.split(' '))[1])
            r[ridx].timevector.append(t)
            r[ridx].addsnapshot();
            sidx = sidx + 1
            r[ridx].snap[sidx].solidcoord = [None]*nsolids
        # Ellipsoid dimensions
        if('#section space' in line):
            spacecoord = [];
            line = fptr.readline()
            if(printstatus):
                print(line)
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                if('ellipse' in line):
                    if(spacecoord):
                        r[ridx].snap[sidx].spacecoord.append(np.array(spacecoord))
                        spacecoord=[];
                    recordStatus = True;
                    if('e' in line[0] and recordStatus):
                        line = line.strip()
                        cstring_tmp = line.split('ellipse')
                        cstring_tmp = cstring_tmp[1].strip()
                        cstring = cstring_tmp.split(' ')
                        spacecoord.append([float(cstring[1]), float(cstring[2]), float(cstring[3])])
                line = fptr.readline()
            # When it exists, if spacecoord has not been recorded, record it.
            if(spacecoord):
                r[ridx].snap[sidx].spacecoord.append(np.array(spacecoord))
                spacecoord=[];
        # Filament coordinates
        if('#section fiber' in line): 
            fcoord=[];
            line = fptr.readline()
            if(printstatus):
                print(line)
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                if('f2' in line or 'f 2' in line):
                    recordStatus = False;
                if('f1' in line or 'f 1' in line):
                    if(fcoord):
                        r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                        fcoord=[];
                    recordStatus = True;
                elif(' ' in line[0] and recordStatus):
                    line = line.strip()
                    cstring = line.split(' ')
                    fcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                line = fptr.readline()
            # When it exists, if fcoord has not been recorded, record it.
            if(fcoord):
                r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                fcoord=[];
        # Crosslinker hands
        if('#section solid' in line and 'filonly' != tag):
            line = fptr.readline()
            line = line.strip()
            if(len(line) == 0):
                line = fptr.readline()
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                # Get solid ID
                line = line.strip()
                cstring = line.split(' ')
                if(len(cstring) == 2):
                    cstring = cstring[0].split(':') # Added for dynamic multimer
                else:
                    cstring = cstring[1].split(':') # changed from [0] to [1]
                solidid = int(cstring[1]) 
                solidcoord = [];
                line = fptr.readline()
                while(line[0]!='d' and not('section' in line)):
                    # Get solid coordinates
                    line = line.strip()
                    cstring = line.split(' ')
                    solidcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                    line = fptr.readline()
                if(solidcoord):
                    r[ridx].snap[sidx].solidcoord[solidid-1] = np.array(solidcoord[0][:])
            # Bound crosslinker hands
            cstring_temp = []
            if('#section single A' in line and 'filonly' != tag):
                line = fptr.readline()
                if(printstatus):
                    print(line)
                while(not('section' in line)):
                    if(printstatus):
                        print(line)
                    if line[0]=='w':
                        line = line.strip()
                        cstring = line.split(' ')
                        if(len(cstring) <= 8):
                            cstring_temp2 = cstring[0].split(':')
                            cstring_temp = float(cstring_temp2[1])
                            if(len(cstring) == 7):
                                r[ridx].snap[sidx].partnerid.append(int(cstring[1][1:len(cstring[1])]))
                                r[ridx].snap[sidx].selfblobid.append(int(cstring[2][1:len(cstring[2])]))
                                r[ridx].snap[sidx].dimerid.append(int(cstring_temp))
                            elif(len(cstring) == 8):
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[1][1:len(cstring[1])]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[3][1:len(cstring[3])]))
                                r[ridx].snap[sidx].boundid.append(int(cstring_temp))
                            else:
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[2][1:len(cstring[2])]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[4][1:len(cstring[4])]))
                        elif(len(cstring) > 8):
                            cstring_temp2 = cstring[1].split(':')
                            cstring_temp = float(cstring_temp2[1])
                            if(len(cstring) == 10):
                                r[ridx].snap[sidx].partnerid.append(int(cstring[3]))
                                r[ridx].snap[sidx].selfblobid.append(int(cstring[5]))
                                r[ridx].snap[sidx].dimerid.append(int(cstring_temp))
                            elif(len(cstring) == 11):
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[3]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[6]))
                                r[ridx].snap[sidx].boundid.append(int(cstring_temp))
                    line = fptr.readline()
    fptr.close()
    print('Number of snapshots='+str(len(r[ridx].snap)))
    fptr.close()
    return r

def getfillength(fc):
    Nbeads = (np.shape(fc))[0]
    L = 0
    for i in range(0,Nbeads-1):
        L = L + np.linalg.norm(fc[i]-fc[i+1])
    return L

def getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename):
    Nsnaps = len(r[0].snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    fvec = np.zeros((2,Ndatapoints))
    fvec[0][:]=np.arange(0,Nsnaps,deltasnap)
    for SREF in range(0,Nsnaps,deltasnap):
        print(SREF,flush=True)
        filcoord = r[0].snap[SREF].filcoord 
        sphericalc = []
        actincounter = ic.generatedensityfield(meshvar, unitnormal, filcoord, False, ic.SearchAlgoType.LISTEDSEARCH, Nlists, sphericalc)
        is_all_zero = np.argwhere((actincounter == 0.0))
        fvec[1][SREF] = 1-len(is_all_zero)/Ntriangles
    pd.DataFrame(fvec, index=['Time','Frac_occupied']).to_csv(outfilename)

def getvaspprops(rset, deltasnap, N, Rval, outputfile):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    bar_min_snap = int(0.95*Nsnaps)
    bar_data_raw = np.array([])
    collectstatus = False
    datawritestatus = False
    datamatrix = np.zeros((19,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    crosslinker_n_fil_bar = [[],[],[],[]]
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        nvasp_tmp = np.array([])
        nvalency_tmp = np.array([])
        Nfil_per_solid_tmp = np.array([])
        if SREF>=bar_min_snap:
            collectstatus = True
        countermat = np.zeros((4,len(rset)))
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                if(len(crossboundblobid)):
                    datawritestatus = True
                unique, counts = np.unique(crossboundblobid, return_counts=True)
                nvasp_tmp = np.append(nvasp_tmp,len(unique))
                nvalency_tmp = np.append(nvalency_tmp,np.array(counts))
                # Get the filament IDs for each solid and calculate Nhands/Nfil
                # This helps you understand the number of hands in each solid that are bound to the same filament.
                # The distribution of Nhands/Nfil - closer to 1 suggests that each bound hand is bound to a different filament
                # >1 suggests that multiple hands of the solid are bound to the same fil ID
                # <1 is not possible.
                for ublobid in enumerate(unique):
                    # Find the locs using argwhere
                    locs = np.argwhere(crossboundblobid==ublobid)
                    # Corresponding filament IDs that share this 
                    # Crosslinker in the bond
                    filid_solid = crossboundfilid[locs]
                    # Get filament IDs and the counts
                    unique_fid, counts_fil = np.unique(filid_solid, return_counts=True)
                    # Get the total number of unique filaments
                    lval = len(unique_fid)
                    countermat[lval-1][runidx] = countermat[lval-1][runidx] + 1
                    Nfil_per_solid_tmp = np.append(Nfil_per_solid_tmp,lval)
                if collectstatus:
                    bar_data_raw = np.append(bar_data_raw, Nfil_per_solid_tmp)
                    for pos in range(0,4):
                        crosslinker_n_fil_bar[pos].append(countermat[pos][runidx])

        meanvec = np.mean(countermat,1)
        sstdvec = np.std(countermat,1)

        if(datawritestatus):
            datamatrix[1][scounter]= np.mean(nvasp_tmp)
            datamatrix[2][scounter]= np.std(nvasp_tmp)
            datamatrix[3][scounter]= np.mean(nvalency_tmp)
            datamatrix[4][scounter]= np.std(nvalency_tmp)
            datamatrix[5][scounter]= np.mean(Nfil_per_solid_tmp)
            datamatrix[6][scounter]= np.std(Nfil_per_solid_tmp)
            datamatrix[9][scounter]= meanvec[0]
            datamatrix[10][scounter]= sstdvec[0]
            datamatrix[11][scounter]= meanvec[1]
            datamatrix[12][scounter]= sstdvec[1]
            datamatrix[13][scounter]= meanvec[2]
            datamatrix[14][scounter]= sstdvec[2]
            datamatrix[15][scounter]= meanvec[3]
            datamatrix[16][scounter]= sstdvec[3]

        scounter = scounter + 1 
        if(len(bar_data_raw)):
            datamatrix[7][0] = np.mean(bar_data_raw)
            datamatrix[8][0] =  np.std(bar_data_raw)
            for pos in range(0,4):
                datamatrix[17][pos] = np.mean(crosslinker_n_fil_bar[pos])
                datamatrix[18][pos] =  np.std(crosslinker_n_fil_bar[pos])
        
    # Write data to file
    print('Saving in file named '+outputfile,flush=True)
    # Nfpc = Number of filaments per crosslinker
    pd.DataFrame(datamatrix, index=['Time', 'Ncross_mean','Ncross_std','Valency_mean','Valency_std',
                                    'Nfil_per_cross_mean','Nfil_per_cross_std',
                                    'Bar_Nfil_per_cross_mean','Bar_Nfil_per_cross_std',
                                    'Mean_Nfpc_1_fil','Std_Npc_1_fil',
                                    'Mean_Nfpc_2_fil','Std_Nfpc_2_fil',
                                    'Mean_Nfpc_3_fil','Std_Nfpc_3_fil',
                                    'Mean_Nfpc_4_fil','Std_Nfpc_4_fil',
                                    'Bar_Mean_Nfpc','Bar_Std_Nfpc'
                                    ]).to_csv(outputfile)

def find_chains_rings_and_doubles(mapping_dict):
    visited = set()
    chains = []
    rings = []
    doubles = set()

    def trace_chain(start_id):
        stack = [(start_id, None)]
        path = []
        ring_detected = False
        ring_start = None
        while stack:
            current_id, prev_id = stack.pop()
            # Rings: Path circled back to a previously recorded link
            if current_id in path:
                ring_start = path.index(current_id)
                ring_detected = True
                break
            # Skip because this path was already recorded
            if current_id in visited:
                continue
            visited.add(current_id)
            path.append(current_id)
            if current_id in mapping_dict:
                next_ids = mapping_dict[current_id]
                # Double bonds: Both hands bound to the same monomer
                if len(next_ids) == 2 and next_ids[0] == next_ids[1]: # Check for double bonds
                    doubles.add(tuple(sorted([current_id, next_ids[0]]))) # Record the double bond
                    return # Does not continue path if double bond detected
                # Continue down the chain
                for next_id in next_ids:
                    if next_id != prev_id: # Avoid going back to the previous node
                        stack.append((next_id, current_id))
        
        if ring_detected:
            rings.append(path[ring_start:] + [path[ring_start]]) # Record the ring
        elif len(path) > 1: 
            chains.append(path) # Record the chain

    for key in mapping_dict.keys():
        if key not in visited:
            trace_chain(key)

    return chains, rings, list(doubles)

def getDynamicMultimerprops(rset, deltasnap, N, Rval, outputfile, outputfile2):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datawritestatus = False
    datamatrix = np.zeros((7,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    length_df = pd.DataFrame()
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The variables below together represent information of each bond in the network.

                # ID of monomer that is bound to a filament
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid

                # ID of monomer that is bound to other monomer(s) --> ID of multimers
                selfblobid = rset[runidx].snap[SREF].selfblobid

                # ID of monomer that this monomer is bound to
                partneridblob = rset[runidx].snap[SREF].partnerid

                # Total number of monomers in the simulation
                nSol = len(rset[runidx].snap[SREF].solidcoord)

                # Count the number of multimers and separate into those bound to 1 or 2 monomers
                multimer_count = Counter(selfblobid)
                Multimer_1 = [num for num, freq in multimer_count.items() if freq == 1]
                Multimer_2 = [num for num, freq in multimer_count.items() if freq == 2]
                
                # Count the number of monomers not bound to other monomers
                all_multimer_blobs = [number for number in range(1,nSol+1)]
                Multimer_0 = [number for number in all_multimer_blobs if number not in selfblobid]

                # Sanity check
                check_total = len(Multimer_0)+len(Multimer_1)+len(Multimer_2)
                if check_total != nSol:
                    print('error: M_0 + M_1 + M_2 != nSolids')

                # Count the number of monomers not bound to a filament
                freeid_blob = [number for number in all_multimer_blobs if number not in crossboundblobid]

                # Sanity check
                check_total = len(freeid_blob)+len(crossboundblobid)
                if check_total != nSol:
                    print('error: bound + free != nSolids')
                
                # Separate all monomers bound to a filament by multimer status
                bound_multimer_0 = [number for number in Multimer_0 if number in crossboundblobid]
                bound_multimer_1 = [number for number in Multimer_1 if number in crossboundblobid]
                bound_multimer_2 = [number for number in Multimer_2 if number in crossboundblobid]

                # Separate all monomers not bound to a filament by multimer status
                free_multimer_0 = [number for number in Multimer_0 if number not in crossboundblobid]
                free_multimer_1 = [number for number in Multimer_1 if number not in crossboundblobid]
                free_multimer_2 = [number for number in Multimer_2 if number not in crossboundblobid]

                # Sanity checks
                check_total_b = len(bound_multimer_0)+len(bound_multimer_1)+len(bound_multimer_2)
                check_total_f = len(free_multimer_0)+len(free_multimer_1)+len(free_multimer_2)
                check_total_t = check_total_b + check_total_f
                if check_total_b != len(crossboundblobid):
                    print('error: sum(bound_multimers) != crossboundblobid')
                if check_total_f != len(freeid_blob):
                    print('error: sum(free_multimers) != freeid_blob')
                if check_total_t != nSol:
                    print('error: sum(multimers) != nSolids')
                else:
                    datawritestatus = True
                
                # Make a dictionary that maps the IDs of each monomer-monomer bond
                partnermap_blob_dict = defaultdict(list)
                for key, value in zip(selfblobid, partneridblob):
                    partnermap_blob_dict[key].append(value)
                partnermap_blob_dict = dict(partnermap_blob_dict)
                
                # Check for different multimer structure types
                chains, rings, doubles = find_chains_rings_and_doubles(partnermap_blob_dict)

                # Get the length of each multimer from each multimer structure type
                chain_lengths = [len(chain) for chain in chains]
                ring_lengths = [len(set(ring)) for ring in rings]
                double_lengths = [len(double) for double in doubles]

                # Sanity check
                for double in doubles:
                    if len(double) != 2:
                        print('Warning: Double has length != 2')
                
                # Make list of all multimer lengths and count the frequency of each multimer length
                all_lengths = chain_lengths + ring_lengths + double_lengths
                length_frequencies = Counter(all_lengths)
                length_df_temp = pd.DataFrame(length_frequencies.items(), columns=['Length', f'{scounter}'], dtype='Int64')
                if length_df.empty:
                    length_df = length_df_temp
                else:
                    length_df = pd.merge(length_df, length_df_temp, on='Length', how='outer')
        
        # Record multimer state
        if(datawritestatus):
            # Does not look across multimer bond
            datamatrix[1][scounter]= len(bound_multimer_0)
            datamatrix[2][scounter]= len(bound_multimer_1)
            datamatrix[3][scounter]= len(bound_multimer_2)
            datamatrix[4][scounter]= len(free_multimer_0)
            datamatrix[5][scounter]= len(free_multimer_1)
            datamatrix[6][scounter]= len(free_multimer_2)

            length_df.fillna(0, inplace=True)
            length_df.set_index('Length', inplace=True)
            length_df.sort_index(inplace=True)

        scounter = scounter + 1 
        
    # Write data to file
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix, index=['Time', 'bound_multimer_0', 'bound_multimer_1', 'bound_multimer_2', 
                                    'free_multimer_0', 'free_multimer_1', 'free_multimer_2']).to_csv(outputfile)
    print('Saving in file named '+outputfile2,flush=True)
    pd.DataFrame(length_df).to_csv(outputfile2)

def getDynamicEllipseprops(rset, deltasnap, N, Rval, outputfile):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datawritestatus = False
    datamatrix = np.zeros((4,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                spacecoord_tmp = np.array(rset[runidx].snap[SREF].spacecoord)[0][0]
                print(spacecoord_tmp)
                # Get ellipsoid dimensions
                if(len(spacecoord_tmp)):
                    datawritestatus = True
                a1 = spacecoord_tmp[0]
                a2 = spacecoord_tmp[1]
                a3 = spacecoord_tmp[2]
        # Record ellipsoid dimensions
        if(datawritestatus):
            datamatrix[1][scounter]= a1
            datamatrix[2][scounter]= a2
            datamatrix[3][scounter]= a3
        scounter = scounter + 1 
    # Write data to file
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix, index=['Time', 'a1', 'a2', 'a3']).to_csv(outputfile)

def getDynamicEllipseFOprops(item_tuple):
    ref, ellipsespan, filcoord = item_tuple
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Snap ="+str(ref)+" Current Time =", current_time, flush=True)
    # Make mesh for this snapshot
    f=4
    meshvar=ic.icosphere(f,ellipsespan*1.025);
    print('icosphere created')
    vertices = meshvar[0]
    triangles = meshvar[1]
    Ntriangles = np.shape(triangles)[0]
    unitnormal = ic.generateunitnormals(meshvar,1.0)
    Nlists = ic.generateNeighborList(meshvar)
    sphericalc = []
    actincounter = ic.generatedensityfield(meshvar, unitnormal, filcoord, True, ic.SearchAlgoType.LISTEDSEARCH, Nlists, sphericalc)
    print('mesh created', flush=True)
    # Generate Fraction Occupied Data
    is_all_zero = np.argwhere((actincounter == 0.0))
    result = 1-len(is_all_zero)/Ntriangles
    return result

def getfillengthprops(rset, deltasnap, N, Rval, outputfile, outputfile2):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    print(Nsnaps,flush=True);
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datamatrix = np.zeros((7,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    scounter = 0
    Lmatrix = []
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        Nfilsnap = []
        L = []
        Lsum = []
        for runidx in range(0,len(rset)):
            Lsum_r = 0
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<nsnaplocal):
                filcoord = (rset[runidx].snap[SREF].filcoord)
                Nfilsnap.append(len(filcoord))
                for f, fc in enumerate(filcoord):
                    Lfil_r = getfillength(fc)
                    Lsum_r = Lsum_r + Lfil_r
                    L.append(Lfil_r)
            Lsum.append(np.sum(Lsum_r))
        Lmatrix.append(L)
        datamatrix[1][scounter]= np.mean(Nfilsnap)
        datamatrix[2][scounter]= np.std(Nfilsnap)
        datamatrix[3][scounter]= np.mean(L)
        datamatrix[4][scounter]= np.std(L)
        datamatrix[5][scounter]= np.mean(Lsum)
        datamatrix[6][scounter]= np.std(Lsum)
        scounter = scounter + 1 
        
    # Write data to file
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix, index=['Time', 'Nfil_mean','Nfil_std','Lfil_mean','Lfil_std',
                                    'Lsum_mean','Lsum_std']).to_csv(outputfile)
    # Calculate moving average
    MOVAVG_NSNAPS = 5
    scounter = 0
    # This matrix will hold the PDF fillength
    maxlen = 6.5 
    nbins = 130 
    fillen_bin_edges = np.linspace(0,maxlen,nbins+1) 
    datamatrix2 = np.zeros((Ndatapoints,nbins)) 
    for SREF in range(0,Nsnaps,deltasnap):
        minsnap = np.max([0,SREF-MOVAVG_NSNAPS])
        maxsnap = np.min([Nsnaps,SREF+MOVAVG_NSNAPS])
        Ltemp = []
        for MSREF in range(minsnap,maxsnap):
            Ltemp = Ltemp + Lmatrix[MSREF]
        kde = gaussian_kde(Ltemp)
        x_plot = fillen_bin_edges[:-1] + (maxlen/nbins)/2
        pdfvec = kde(x_plot)
        datamatrix2[scounter][:] = pdfvec
        scounter = scounter + 1
    pd.DataFrame(datamatrix2, index=np.arange(0,Nsnaps,deltasnap)).to_csv(outputfile2)

def geteigensnap(filcoord):
    snapcoordmat=np.array([[]])
    for f, fc in enumerate(filcoord):
        Nbeads = (np.shape(fc))[0]
        filcoordmat = np.zeros((Nbeads,3))
        for i in range(0,Nbeads):
            filcoordmat[i,:] = np.array(fc[i,:])
        if f==0:
            snapcoordmat = filcoordmat.copy()
        else:
            axisval = 0
            snapcoordmat = np.concatenate((snapcoordmat,filcoordmat),axis=axisval)

    COM = np.mean(snapcoordmat,0)
    snapcoordmat = (snapcoordmat - COM)
    crosscorr  = np.matmul(np.transpose(snapcoordmat),snapcoordmat)/len(snapcoordmat)
    eval, evec =np.linalg.eig(crosscorr)
    eval = np.flip(np.sort(eval))
    eval = 2*np.sqrt(eval)
    return eval

def geteigentrajectory(rset, deltasnap, N, Rval, outputfile):
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    # Row represents data 
    datamatrix = np.zeros((4,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        PC1 = []; PC2 = []; PC3 = [];
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                filcoord = (rset[runidx].snap[SREF].filcoord)
                PCeigen = geteigensnap (filcoord)
                PC1.append(PCeigen[0])
                PC2.append(PCeigen[1])
                PC3.append(PCeigen[2])
        datamatrix[1][scounter] = np.mean(PC1)
        datamatrix[2][scounter] = np.mean(PC2)
        datamatrix[3][scounter] = np.mean(PC3)
        scounter = scounter + 1
    # Write data to file
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix,['Time','eig1','eig2','eig3']).to_csv(outputfile)

def analyzetrajectory(fpathvar,N, dirname, outfilename):
    f=4
    print(dirname,flush=True)
    meshvar=ic.icosphere(f,1.0);
    vertices = meshvar[0]
    triangles = meshvar[1]
    Ntriangles = np.shape(triangles)[0]
    unitnormal = ic.generateunitnormals(meshvar,1.0)
    Nlists = ic.generateNeighborList(meshvar)
    print('mesh created', flush=True)
    deltasnap = 1
    dirlist = []
    dirlist.append(dirname)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename)

######## Frontend Functions ########

def frontend_FO_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/FO_'+foldername+'_'+dirname+'.csv'
    print(foldername, flush=True)
    analyzetrajectory(fpathvar,N, dirname, outputfile)
    print("The End....", flush=True)

def frontend_cross_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getvaspprops(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_cross_DynamicMultimer_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    outputfile2 = args[1]+'outputfiles/Length_Freq_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getDynamicMultimerprops(rset, 1, N, Rval, outputfile, outputfile2)
    print("The End....", flush=True)

def frontend_dynamic_ellipse_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Dynamic_Ellipse_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getDynamicEllipseprops(rset, 1, N, Rval, outputfile)
    print("The End....", flush=True)

def frontend_dynamic_ellipse_FO_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/FO_'+foldername+'_'+dirname+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    read = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(read[0])
    print('Read repid '+str(repid),flush=True)
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    deltasnap = N + 1
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    fvec = np.zeros((2,Ndatapoints))
    fvec[0][:] = np.arange(0, Nsnaps, deltasnap)
    # Make tuple list of ellipsoid dimensions and filament coordinates
    tuple_list = []
    for ref in range(0,Nsnaps,deltasnap):
        ellipsespan = np.array(r.snap[ref].spacecoord)[0][0]
        filcoord = r.snap[ref].filcoord
        tuple_list.append((ref, ellipsespan, filcoord))
    # Set up parallel jobs
    res = Parallel(n_jobs=30)(delayed(getDynamicEllipseFOprops)(item_tuple) for item_tuple in tuple_list)
    fvec[1][:] = np.array(res)
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(fvec, index=['Time','Frac_occupied']).to_csv(outputfile)
    print("The End....", flush=True)

def frontend_filprops_Rval(N, Rval, Nreps, *args):
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Fillength_'+foldername+'R_'+str(Rval)+'.csv'
    outputfile2 = args[1]+'outputfiles/Fillength_dist_'+foldername+'R_'+str(Rval)+'.csv'
    print(foldername, flush=True)
    rset = []
    for repid in range(0,Nreps):
        dirname = 'R_'+str(Rval)+'_r_'+str(repid)
        print('Reading trajectory '+dirname,flush=True)
        r = readcytosimtraj(fpathvar+dirname+'/','filonly')
        print('Trajectory read..',flush=True)
        rset.append(r[0])
        print('Read repid '+str(repid),flush=True)
    getfillengthprops(rset, 1, N, Rval,outputfile, outputfile2)
    print("The End....", flush=True)

def frontend_eigval_Rval_repid(N, Rval, repid, *args):
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Eigval'+foldername+'R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)
    rset = []
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    geteigentrajectory(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)