# Eliminate only loops longer than X

from Bio import PDB
from Bio.PDB.DSSP import *
from sys import exit
import numpy as np
import argparse

# DSSP codes
# H = alpha-helix
# B = residue in isolated beta-bridge
# E = extended strand, participates in beta ladder
# G = 3-helix (310 helix)
# I = 5 helix (pi-helix)
# T = hydrogen bonded turn
# S = bend
ss_dic = {'H':0, 'B':1, 'E':2, 'G':3, 'I':4, 'T': 5, 'S':6, '-':7}

def get_resi_pdb():
    # Get the residue numbers in pdb file
    p  = PDB.PDBParser()
    s = p.get_structure(name, inputs.pdb)

    rr = s.get_residues()
    resi = {}
    for r in rr:
        n1, n2, n3 = r.get_id()
        ch = r.get_parent().get_id()
        m = r.resname
        print(n2, m)
        if ch in resi.keys():
            resi[ch].append(n2)
        else:
            resi[ch] = [n2]

    for k,v in resi.items():
        resi[k] = np.sort(v)
    return resi
    
################################################
parser = argparse.ArgumentParser(description='Remove loops form a pdb. Options -pdb, -loop_length, -h')
parser.add_argument('-pdb', action="store", dest="pdb", help="Input pdb")
parser.add_argument('-loop_length', action="store", dest="loop_length", help="Minimum number of residues un loop")

inputs = parser.parse_args()
if inputs.pdb == None:
    print('Usage: remove_loops.py -pdb pdb_in -loop_length')
    exit()
if inputs.loop_length == None:
    inputs.loop_length = 12
################################################

# Check if given pdb id file or complete path
if len(inputs.pdb.split('/')) > 0:
    pdb_name = inputs.pdb.split('/')[-1]
    pdb_path = inputs.pdb[0:-len(pdb_name)] 
else:
    pdb_name = inputs.pdb
    pdb_path = './'
name = pdb_name.split('.')[0]



# Get residue numbers
seqid = get_resi_pdb()

# Compute secondary structure using DSSP
dssp_dict=dssp_dict_from_pdb_file(inputs.pdb,DSSP='/Users/iecheverria/SOFTW/dssp-2.0.4/mkdssp')

# Go over residues to get SS
#ss = np.empty([len(seqid),3],dtype=int)

ss_vals = {}
for k, v in seqid.items():
    print(k)
    for r in v:
        ss_pred = dssp_dict[0][(k, (' ', r, ' '))][1]
        val = ss_dic[ss_pred]
        if k not in ss_vals.keys():
            ss_vals[k] = [[r, val, 0]]
        else:
            ss_vals[k].append([r, val, 0])

for k, resi in ss_vals.items():
    ss_vals[k] = np.array(resi)
    

# Find loop longer than loop_length
for k, resi in ss_vals.items():
    i = 0
    while i < len(resi[:,0])-1:
        if resi[i,1] in [0,1,2,3,4]:
            resi[i,2] = 0
            i += 1
        else:
            m = 0
            ss_run = resi[i,1]
            while ss_run not in [0,1,2,3,4] and (i + m) < len(resi[:,0])-1:
                m += 1
                ss_run = resi[(i+m),1]
            print(m, int(inputs.loop_length))
            if m >= int(inputs.loop_length):
                resi[i:(i+m),2] = inputs.loop_length
            else:
                resi[i:(i+m),2] = resi[i:(i+m),1]
            i = i + m

print(ss_vals['A'][:,0])           

# Write pdb
out = open(name+'.nl.pdb','w')
for line in open(inputs.pdb):
    vals = line.split()
    if vals[0]== 'ATOM':
        chain = line[21:22]
        resi = int(line[22:26])
        seq_struct = ss_vals[chain][ss_vals[chain][:,0]==resi,:]
        if seq_struct[0][2] != int(inputs.loop_length):
            out.write(line)
        
    else:
        out.write(line)

out.close()
