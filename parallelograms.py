from Bio.PDB import *
import math
import csv
import pandas as pd

def distance(atom1,atom2):
    return math.sqrt(math.pow(atom1[0]-atom2[0],2)+math.pow(atom1[1]-atom2[1],2)+math.pow(atom1[2]-atom2[2],2))

if __name__ == '__main__':
    length  = 161
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
        flag = 0
        for atom in generator:
            if atom.get_name() == "CA":
                CA_list.append(atom.get_coord())
            if atom.get_name() == "CB":
                CB_list.append(atom.get_coord())
                flag = 1
        #  handle glycine
        if flag == 0:
            print(i)
            CB_list.append("GLY")

    print(len(CA_list), len(CB_list))
    w, h = length, length
    ACAD = [[0 for x in range(w)] for y in range(h)]
    ACBC = [[0 for x in range(w)] for y in range(h)]
    ACBD = [[0 for x in range(w)] for y in range(h)]
    ADBC = [[0 for x in range(w)] for y in range(h)]
    ADBD = [[0 for x in range(w)] for y in range(h)]
    CDBD = [[0 for x in range(w)] for y in range(h)]
    for i in range(length):
        for j in range(length):
            if CB_list[i] == "GLY" or CB_list[j] == "GLY":
                #print(i,j)
                ACAD[i][j] = "GLY"
                ACBC[i][j] = "GLY"
                ACBD[i][j] = "GLY"
                ADBC[i][j] = "GLY"
                ADBD[i][j] = "GLY"
                CDBD[i][j] = "GLY"
            elif i != j:
                p = list()
                a = CA_list[i]
                b = Vector(CB_list[i] - a)
                A = Vector(a)
                B = A + b/Vector.norm(b)
                c = CA_list[j]
                d = Vector(CB_list[j] - c)
                C = Vector(c)
                D = C + d/Vector.norm(d)
                AB = Vector.norm(B - A)
                AC = Vector.norm(C - A)
                AD = Vector.norm(D - A)
                BC = Vector.norm(C - B)
                BD = Vector.norm(D - B)
                CD = Vector.norm(D - C)
                ACAD[i][j] = AC/AD
                ACBC[i][j] = AC/BC
                ACBD[i][j] = AC/BD
                ADBC[i][j] = AD/BC
                ADBD[i][j] = AD/BD
                CDBD[i][j] = BC/BD
                #(AB, AC, AD, BC, BD, and CD)
                #(AC/AD, AC/BC, AC/BD, AD/BC, AD/BD, BC/BD)

    with open("ACADout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(ACAD)

    with open("ACBCout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(ACBC)

    with open("ACBDout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(ACBD)

    with open("ADBCout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(ADBC)

    with open("ADBDout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(ADBD)

    with open("CDBDout.csv","w") as f:
        wr = csv.writer(f,delimiter="\n")
        wr.writerow(CDBD)
