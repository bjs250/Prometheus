from Bio.PDB import *
import math
import csv
import pandas as pd

def distance(atom1,atom2):
    return math.sqrt(math.pow(atom1[0]-atom2[0],2)+math.pow(atom1[1]-atom2[1],2)+math.pow(atom1[2]-atom2[2],2))

if __name__ == '__main__':
    length  = 10
    parser = PDBParser()
    structure = parser.get_structure('protein', '3fa0.pdb')
    residues = list(structure.get_residues())

    # create list of residue atoms and their coordinates
    res_atoms = list()
    res_coords = list()
    for i in range(length):
        residue = residues[i]
        atoms = list()
        coords = list()
        for atom in residue:
            atoms.append(atom)
            coords.append(atom.get_coord())
        res_atoms.append(atoms)
        res_coords.append(coords)

    # create list of CA and CB coordinates
    CA_list = list()
    CB_list = list()
    for i in range(length):
        generator = residues[i].get_atoms()
        for atom in generator:
            if atom.get_name() == "CA":
                CA_list.append(atom.get_vector())
            if atom.get_name() == "CB":
                CB_list.append(atom.get_vector())

    w, h = length, length
    AngleMatrix = [[0 for x in range(w)] for y in range(h)]
    for i in range(length):
        for j in range(length):
            if i != j:
                AngleMatrix[i][j] = calc_dihedral(CA_list[i], CB_list[i], CB_list[j], CA_list[j]) * 180 / math.pi

    Matrix = [[0 for x in range(w)] for y in range(h)]
    for i in range(length):
        for j in range(length):
            if i != j:
                min_dist = list()
                for atom1 in res_atoms[i]:
                    dist = list()
                    for atom2 in res_atoms[j]:
                        dist.append(distance(atom1.get_coord(),atom2.get_coord()))
                    min_dist.append(min(dist))
                Matrix[i][j] = min(min_dist)


    with open("out.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(Matrix)

    with open("aout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(AngleMatrix)
