#!/usr/bin/python

import os
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        help = "input PDB file",
                        required = True)
    parser.add_argument('-c',
                        help = "chain to include in the distance calculation",
                        required = True)
    parser.add_argument('-lim',
                        help = "maximum cutoff for distances (in angstorms)",
                        default = "5")
    parser.add_argument('-atm',
                        help = 'atom type for which distances are to be calculated: C or H',
                        required = True)
    args = parser.parse_args()

    atom_lines = []
    with open('%s/%s' %(os.getcwd(), args.i)) as input_file:
        for line in input_file:
            line = line.strip()
            line = line.split()
            if ((line[0] == 'ATOM') and (line[4] == args.c) and (line[11] == args.atm)):
                atom_lines.append(line)
    distances = []
    for atom_line in atom_lines:
        coordinates_1 = np.array((float(atom_line[6]), float(atom_line[7]), float(atom_line[8])))
        for j in range(len(atom_lines)):
            coordinates_2 = np.array((float(atom_lines[j][6]), float(atom_lines[j][7]), float(atom_lines[j][8])))
            distance = np.sqrt(np.sum((coordinates_1 - coordinates_2)**2))
            distance = np.round(distance, decimals=3)
            if ((distance <= float(args.lim)) and (distance > 0.0)):
                distances.append([atom_line[5], atom_line[2], atom_lines[j][5], atom_lines[j][2], distance])
        atom_lines.remove(atom_line)
    for distance in distances:
        print "%s,%s,%s,%s,%s" %(distance[0], distance[1], distance[2], distance[3], distance[4])

if __name__ == '__main__':
    main()
