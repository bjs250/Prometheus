from Bio.PDB import *
import math
import csv
import pandas as pd

if __name__ == '__main__':
    length  = 50
    parser = PDBParser()
    structure = parser.get_structure('protein', '3fa0.pdb')
    residues = list(structure.get_residues())
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

    CA_list = list()
    CB_list = list()

    for i in range(length):
        generator = residues[i].get_atoms()
        for atom in generator:
            if atom.get_name() == "CA":
                CA_list.append(atom.get_vector())
            if atom.get_name() == "CB":
                CB_list.append(atom.get_vector())

    print(calc_dihedral(CA_list[2],CB_list[2],CB_list[6], CA_list[6]) * 180 / math.pi)
    print(calc_dihedral(CA_list[3],CB_list[3],CB_list[7], CA_list[7]) * 180 / math.pi)
