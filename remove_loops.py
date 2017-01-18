#!/usr/bin/env python

# Eliminate loops longer than loop_length from a pdb
# iecheverria - Sali lab - UCSF
# Uses Biopython, DSSP and argparse

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
        m = r.resname
        resi[n2] = m

    seqid = np.sort(np.array(resi.keys()))
    return seqid
    
################################################
parser = argparse.ArgumentParser(description='Remove loops form a pdb. Options -pdb, -loop_length, -h')
parser.add_argument('-pdb', action="store", dest="pdb", help="Input pdb")
parser.add_argument('-loop_length', action="store", dest="loop_length", help="Minimum number of residues un loop")

inputs = parser.parse_args()
if inputs.pdb == None:
    print(' ')
    print('Usage: remove_loops.py -pdb pdb_in -loop_length')
    print(' ')
    exit()
if inputs.loop_length == None:
    inputs.loop_length = 8
################################################

# Check if inputs.pdb is a file or path
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
dssp_dict=dssp_dict_from_pdb_file(inputs.pdb,DSSP='mkdssp')

# Go over residues to get SS
ss = np.empty([len(seqid),3],dtype=int)
for i in range(len(seqid)):
    ss[i,0] = seqid[i]
    ss[i,1] = ss_dic[dssp_dict[0][('A', (' ', seqid[i], ' '))][1]]

# Find loop longer than loop_length
i = 0
while i < len(seqid)-1:
    if ss[i,1] in [0,1,2]:
        ss[i,2] = 0
        i += 1
    else:
        k = 0
        ss_run = ss[i,1]
        while ss_run not in [0,1,2] and (i + k) < len(seqid)-1:
            k += 1
            ss_run = ss[(i+k),1]
        if k >=8:
            ss[i:(i+k),2] = inputs.loop_length
        else:
            ss[i:(i+k),2] = ss[i:(i+k),1]
        i = i + k

# Write pdb
out = open(name+'.nl.pdb','w')
for line in open(inputs.pdb):
    vals = line.split()
    if vals[0]== 'ATOM':
        resi = int(line[22:26])
        seq_struct = ss[ss[:,0]==resi,:]
        if seq_struct[0][2] != 8:
            out.write(line)

out.close()
