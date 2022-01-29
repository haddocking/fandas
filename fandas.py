#!/usr/bin/python

#AUTHOR: Siddarth Narasimhan, NMR Spectroscopy Group, Utrecht University

#all imports
import os
import sys
import argparse
import numpy as np
import shutil

#global definitions for fandas
w_dir = os.getcwd()

amino_acids = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I',
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

secondary_structures = ['a', 'b', 'c', 'n'] #a=alpha helix, b=beta sheet, c=random coil, n=bmrb collected average

atoms = ['H', 'HA', 'HA2', 'HA3', 'HB', 'HB2', 'HB3', 'HG1', 'HG2', 'HG21',
         'HG22', 'HG3', 'HD', 'HD1', 'HD2', 'HD21', 'HD22', 'HD3', 'HE',
         'HE1', 'HE12', 'HE13', 'HE2', 'HE3', 'HH', 'HH11', 'HH12', 'HH2',
         'HH21', 'HH22', 'HZ', 'HZ2', 'HZ3', 'N', 'NG1', 'NG2', 'ND', 'ND1',
         'ND2', 'NH1', 'NH2', 'NZ', 'C', 'CA', 'CB', 'CG', 'CG1', 'CG2', 'CD',
         'CD1', 'CD2', 'CD3', 'CE', 'CE1', 'CE2', 'CH2', 'CZ', 'CZ2', 'CZ3']

atom_type = ['H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H',
             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'N', 'N', 'N', 'N', 'N', 'N',
             'N', 'N', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C',
             'C', 'C', 'C', 'C', 'C', 'C', 'C']

atom_positions  = ['N', 'A', 'A2', 'A3', 'B', 'B2', 'B3', 'G1', 'G2', 'G2',
                   'G2', 'G3', 'D', 'D1', 'D2', 'D2', 'D2', 'D3', 'E',
                   'E1', 'E1', 'E1', 'E2', 'E3', 'H', 'H1', 'H1', 'H2',
                   'H2', 'H2', 'Z', 'Z2', 'Z3', 'N', 'G1', 'G2', 'D', 'D1',
                   'D2', 'H1', 'H2', 'Z', 'C', 'A', 'B', 'G', 'G1', 'G2', 'D',
                   'D1', 'D2', 'D3', 'E', 'E1', 'E2', 'H2', 'Z', 'Z2', 'Z3']

simple_atom_positions  = ['N', 'A', 'A', 'A', 'B', 'B', 'B', 'G', 'G', 'G',
                          'G', 'G', 'D', 'D', 'D', 'D', 'D', 'D', 'E',
                          'E', 'E', 'E', 'E', 'E', 'H', 'H', 'H', 'H',
                          'H', 'H', 'Z', 'Z', 'Z', 'N', 'G', 'G', 'D', 'D',
                          'D', 'H', 'H', 'Z', 'C', 'A', 'B', 'G', 'G', 'G', 'D',
                          'D', 'D', 'D', 'E', 'E', 'E', 'H', 'Z', 'Z', 'Z']

h_atm_ind = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]

n_atm_ind = [33, 34, 35, 36, 37, 38, 39, 40, 41]
    
c_atm_ind = [42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58]

list_2d = ["NH", "HN", "CH", "HC", "HH", "DQSQ", "CC_SPINDIFF_INTRA", "CC_SPINDIFF_INTER", "NCA", 
           "NCO", "NCACX", "NCACX_INTER", "NCOCX", "NCOCA_CB", "CANH", "CONH", "CACONH", "COCANH", 
           "NCAH"]

list_3d = ["NCACX", "NCACX_INTER", "NCOCX", "NCOCA_CB", "SQSQSQ_INTER", "DQSQSQ_INTRA",
           "DQSQSQ_INTER", "CANH", "CONH", "CACONH", "COCANH", "NCAH"]

list_2dd = ["CC_SPINDIFF","HH", "CHHC", "NHHC", "CHH", "NHH", "HHC", "NCACX", "NCOCX"]

list_3dd = ["CHH", "NHH", "NCACX", "NCOCX"]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        help = "input sequence as a text file (REQUIRED)",
                        required = True)
    parser.add_argument('-ss',
                        help = "secondary structures as a text file: 'a'- alpha helix, 'b'- beta sheet & 'c'- random coil; if unspecified, will use BMRB averages: 'n'",
                        default = '')
    parser.add_argument('-bt',
                        help = "BMRB tables as a text file in NMR Star format",
                        default = '')
    parser.add_argument('-btc',
                        help = 'column numbers for residue number, atom name and chemical shifts in the BMRB tables',
                        nargs = 3,
                        default = [1, 3, 4])
    parser.add_argument('-ls',
                        help = 'labelling scheme. Default is Fully Labelled. Other options: fw = Forward Labelled, rv = Reverse Labelled, gl13 = 1,3 Glycerol Labelling, gl2= 2 Glycerol Labelling')
    parser.add_argument('-dl',
                         help = 'list of 13C & 15N (for forward labelling) or 12C & 14N (for reverse labelling) amino acids',
                         nargs = '+',
                         default = '')
    parser.add_argument('-cl',
                         help = 'list of 13C (for forward labelling) or 12C (for reverse labelling) amino acids',
                         nargs = '+',
                         default = '')
    parser.add_argument('-nl',
                         help = 'list of 15N (for forward labelling) or 14N (for reverse labelling) amino acids',
                         nargs = '+',
                         default = '')
    parser.add_argument('-fd',
                         help = 'fractionally deuterated, if you use this flag, it would be implemented automatically',
                         action = 'store_true')
    parser.add_argument('-exp_2d',
                         help = 'list of 2D experiments: %s' %list_2d,
                         nargs = '+',
                         default = '')
    parser.add_argument('-exp_2dd',
                         help = 'list of distance encoded (distance list, -dl and limit -dlm  must be provided) 2D experiments: %s' %list_2dd,
                         nargs = '+',
                         default = '')
    parser.add_argument('-exp_3d',
                         help = 'list of 3D experiments: %s' %list_3d,
                         nargs = '+',
                         default = '')
    parser.add_argument('-exp_3dd',
                         help = 'list of distance encoded (distance list, -dl and limit -dlm must be provided) 3D experiments: %s' %list_3dd,
                         nargs = '+',
                         default = '')
    parser.add_argument('-sl',
                        help = 'automatically assign peak labels in the sparky file',
                        action = 'store_true')
    parser.add_argument('-dlist',
                         help = 'distance list (as a file) in angstorms as follows: "resi_num,atm_nam,resi_num_2,atm_nam_2,dist". EXAMPLE:"2,CA,4,CB,4.5"',
                         default = '')
    parser.add_argument('-dlim',
                         help = 'distance limit in angstorms, default: 5 Angstorm',
                         default = 5)
    parser.add_argument('-o',
                         help = 'names for output, default: "fandas_output"',
                         default = "fandas_output")
    args = parser.parse_args()

    ###(1) give sensible names to form inputs and parse inputs
    global project_name
    global sequence
      ###STEP 1: READ INPUT SEQUENCE AND CHECK THE SEQUENCE USING check_user_input FUNCTION###
    project_name = args.o
    with open("%s/%s" %(w_dir, args.i), 'r') as input_file:
        sequence = list(input_file.readline()) #Parse sequence
    
    for residue in sequence:
        if residue == '\n' : sequence.remove('\n') #Remove new line charecter
        elif residue == ' ' : sequence.remove(' ') #Remove spaces
    
    for i in range(len(sequence)):
        sequence[i] = sequence[i].upper()
    check_user_input(sequence, amino_acids, "Sequence")

    ###STEP 2: ASSIGN SINGLE LETTER SECONDARY STRUCTURE FOR THE RESIDUES###
    ##If the secondary structure is not given, BMRB average shifts (as on 21/07/16) will be assigned##
    if args.ss == '':
        sec_struc = ['n']*len(sequence)
    else:
        with open("%s/%s" %(w_dir, args.ss), 'r') as ss_file:
            sec_struc = list(ss_file.readline())
        for ss in sec_struc:
            if ss == '\n': sec_struc.remove('\n')
            if ss == ' ': sec_struc.remove(' ')
        for i in range(len(sec_struc)):
            sec_struc[i] = sec_struc[i].lower()
        ##Assign bmrb average to make up for missing ss. If the ss assi
        if len(sec_struc) < len(sequence):
            difference = len(sequence) - len(sec_struc)
            for i in range(difference):
                sec_struc.append('n')
        ##If the user defined secondary structure is longer than sequence, chop it##
        elif len(sec_struc) > len(sequence):
            sec_struc = sec_struc[:len(sequence)]
        check_user_input(sec_struc, secondary_structures, "Secondary Structure")
     
    ###(4) assign chemical shifts
    ##assign average shifts ased on the above single letter assignments##
    chem_shifts = assign_chemical_shifts(sequence, sec_struc)

    ###(5) replace the average shifts with provided shifts in the form of BMRB table
    ##Replace the average shifts with provided shifts in the form of BMRB table##
    if args.bt != '':
        chem_shifts = replace_bmrb(chem_shifts, args.bt, args.btc)

    ###(6) incorporate forward, reverse & glycerol labelling schemes
    #configure labelled amino acids list if fwd or rv labelling schemes are used
    labelling_scheme = args.ls
    dl = args.dl
    cl = args.cl
    nl = args.nl
    if ((labelling_scheme == 'fw') or (labelling_scheme == 'rv')):
        fdl = []
        fcl = []
        fnl = []
        rdl = []
        rcl = []
        rnl = []
        if (dl != []):
            if (labelling_scheme == 'fw'):
                fdl = dl
            else:
                rdl = dl
        if (cl != []):
            if (labelling_scheme == 'fw'):
                fcl = cl
            else:
                rcl = cl
        if (nl != []):
            if (labelling_scheme == 'fw'):
                fnl = nl
            else:
                rnl = nl
        #clone amino_acid list onto aa_rm from which all the fwd_labelled amino acids will be removed
        if (labelling_scheme == 'fw'):
            aa_rm = amino_acids
        #remove all the fwd_cn labelled amino acids from the aa_rm list
        if (fdl != []):
            for fw_double in fdl:
                aa_rm.remove(fw_double.upper())
        #remove all the 13C shifts for fwd_n labelled
        #remove all the 15N labelled amino acids from the aa_rm list
        if (fnl != []):
            for aa in fnl:
                chem_shifts = rev_label(chem_shifts, 'rev_c', aa.upper())
                aa_rm.remove(aa.upper())
        #remove all the 15N shifts for fwd_c labelled
        #remove all the 13C labelled amino acids from the aa_rm list
        if (fcl != []):
            for aa in fcl:
                chem_shifts = rev_label(chem_shifts, 'rev_n', aa.upper())
                aa_rm.remove(aa.upper())
        #remove all the 13C & 15N shifts for rev_cn labelled
        if (rdl != []):
            for aa in rdl:
                chem_shifts = rev_label(chem_shifts, 'rev_cn', aa.upper())
        #remove all the 15N shifts for rv_n laabelled
        if (rnl != []):
            for aa in rnl:
                chem_shifts = rev_label(chem_shifts, 'rev_n', aa.upper())
        #remove all the 13C shifts for fwd_c labelled
        if (rcl != []):
            for aa in rcl:
                chem_shifts = rev_label(chem_shifts, 'rev_c', aa.upper())
        #remove all the 13C & 15N shifts for everything in the aa_rm list
        if ((fdl != []) or (fcl != []) or (fnl != [])):
            for aa in aa_rm:
                chem_shifts = rev_label(chem_shifts, 'rev_cn', aa.upper())
    elif labelling_scheme == 'gl13':
        chem_shifts = glycerol_label(sequence, chem_shifts, 13)
    elif labelling_scheme == 'gl2':
        chem_shifts = glycerol_label(sequence, chem_shifts, 2)

    ###7 TAKE FRACTIONAL DEUTERATION INTO ACCOUNT###
    if args.fd == True:
        chem_shifts = fractional_deuteration(sequence, chem_shifts)
    chem_shifts = np.around(chem_shifts, decimals = 4)
    
    ###(8) Read distance lists
    if (args.dlist != ''):
        dlim = args.dlim
        distances = []
        dist_list = []        
        with open("%s/%s" %(w_dir, args.dlist), 'r+') as list_file:
            for line in list_file:
                line = line.rstrip()
                line = line.replace(" ", "")
                line = line.split(',')
                if float(line[4]) <= float(dlim): 
                    distances.append(line[4])
                    dist_list.append([int(line[0]), line[1], int(line[2]), line[3]])
    else:
        dist_list = []
    
    ###(9) Round off the final chemical shifts to 2 decimal places and include options for sparky labels
    chem_shifts = np.around(chem_shifts, decimals = 2)
    global sl
    sl = args.sl

    ###(10) Make Predictions
    exp_2d = args.exp_2d
    exp_2dd = args.exp_2dd
    exp_3d = args.exp_3d
    exp_3dd = args.exp_3dd
    if exp_2d != []:
        for experiment in exp_2d:
            experiment = experiment.upper()
            if experiment == "NH": peaks_proton_heavy(chem_shifts, 'H', 'N', 'H')
            elif experiment == "HN": peaks_proton_heavy(chem_shifts, 'H', 'N', 'N')
            elif experiment == "CH": peaks_proton_heavy(chem_shifts, 'H', 'C', 'C')
            elif experiment == "HC": peaks_proton_heavy(chem_shifts, 'H', 'C', 'H')
            elif experiment == "HH": hh(chem_shifts) 
            elif experiment == "DQSQ": dqsq(chem_shifts)
            elif experiment == "CC_SPINDIFF_INTRA": cc_spindiff_intra(chem_shifts)
            elif experiment == "CC_SPINDIFF_INTER": cc_spindiff_inter(chem_shifts)
            elif experiment == "NCA": nca(chem_shifts)
            elif experiment == "NCO": nco(chem_shifts)
            elif experiment == "NCACX": ncacx(chem_shifts, 2)
            elif experiment == "NCACX_INTER": ncacx_inter(chem_shifts, 2)
            elif experiment == "NCOCX": ncocx(chem_shifts, 2)
            elif experiment == "NCOCA_CB": ncoca_cb(chem_shifts, 2)
            elif experiment == "CANH": canh(chem_shifts, 2)
            elif experiment == "CONH": conh(chem_shifts, 2)
            elif experiment == "CACONH": caconh(chem_shifts, 2)
            elif experiment == "COCANH": cocanh(chem_shifts, 2)
            elif experiment == "NCAH": ncah(chem_shifts, 2)
            else: print "%s is not a valid experiment" %experiment
    
    if exp_3d != []:
        for experiment in exp_3d:
            experiment = experiment.upper()
            if experiment == "NCACX": ncacx(chem_shifts, 3)
            elif experiment == "NCACX_INTER": ncacx_inter(chem_shifts, 3)
            elif experiment == "NCOCX": ncocx(chem_shifts, 3)
            elif experiment == "NCOCA_CB": ncoca_cb(chem_shifts, 3)
            elif experiment == "SQSQSQ_INTER": sqsqsq_inter(chem_shifts)
            elif experiment == "DQSQSQ_INTRA": dqsqsq_intra(chem_shifts)
            elif experiment == "DQSQSQ_INTER": dqsqsq_inter(chem_shifts)
            elif experiment == "CANH": canh(chem_shifts, 3)
            elif experiment == "CONH": conh(chem_shifts, 3)
            elif experiment == "CACONH": caconh(chem_shifts, 3)
            elif experiment == "COCANH": conh(chem_shifts, 3)
            elif experiment == "NCAH": ncah(chem_shifts, 3)
            else: print "%s is not a valid experiment" %experiment

    if ((exp_2dd != []) and (dist_list != [])):
        for experiment in exp_2dd:
            experiment = experiment.upper()
            if experiment == "CC_SPINDIFF": cc_spindiff_dist(chem_shifts, dist_list, distances)
            elif experiment == "HH": hh_dist(chem_shifts, dist_list, distances)
            elif experiment == "CHHC": chhc_dist(chem_shifts, dist_list, distances)
            elif experiment == "NHHC": nhhc_dist(chem_shifts, dist_list, distances)
            elif experiment == "CHH": chh_dist(chem_shifts, dist_list, distances)
            elif experiment == "NHH": nhh_dist(chem_shifts, dist_list, distances)
            elif experiment == "HHC": hhc_dist(chem_shifts, dist_list, distances)
            elif experiment == "NCACX": ncacx_dist(chem_shifts, 2, dist_list, distances)
            elif experiment == "NCOCX": ncocx_dist(chem_shifts, 2, dist_list, distances)
            else: print "%s is not a valid experiment" %experiment

    if ((exp_3dd != []) and (dist_list != [])):
        for experiment in exp_3dd:
            if experiment == "CHH": chh_dist_3d(chem_shifts, dist_list, distances)
            elif experiment == "NHH": nhh_dist_3d(chem_shifts, dist_list, distances)
            elif experiment == "NCACX": ncacx_dist(chem_shifts, 3, dist_list, distances)
            elif experiment == "NCOCX": ncocx_dist(chem_shifts, 3, dist_list, distances)
            else: print "%s is not a valid experiment" %experiment

def peaks_proton_heavy(chem_shifts, atm_1, atm_2, direct):
    #this function produces peak list for CH, HC, HN and NH type experiments
    #atm_2 is by default the heavy atom
    #direct defines the atom on the direct dimension
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atom in enumerate(atom_type): #loop through all atomtypes
            if atom == atm_2: #check if it is a heavy atom
                for h_ind in h_atm_ind: #loop on proton indices
                    if atom_positions[h_ind] == atom_positions[j] and chem_shifts[i][h_ind] != 0 and chem_shifts[i][j] != 0:
                            shift_1.append([residue, i+1, chem_shifts[i][h_ind], atoms[h_ind]])
                            shift_2.append([residue, i+1, chem_shifts[i][j], atoms[j]])
    if atm_1 == direct:
        extension = "%s%s" %(atm_1, atm_2)
        write_2d(shift_1, shift_2, extension.lower())
    else: 
        extension = "%s%s" %(atm_2, atm_1)
        write_2d(shift_2, shift_1, extension.lower())

def hh(chem_shifts):
    #this function produces peak list for intra-residue HH experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atom in enumerate(atom_type): #loop through all atomtypes
            if atom == 'H': #check if it is a carbon atom
                for h_ind in h_atm_ind: #loop through all the carbon indices
                    if h_ind != j: #ensure that the two indices "j" and "h_ind" are not the same
                        if chem_shifts[i][h_ind] != 0 and chem_shifts[i][j] != 0: #check if the chemical shifts of the two are not zero
                            shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                            shift_2.append([residue, i+1, chem_shifts[i][h_ind], atoms[h_ind]])
    write_2d(shift_1, shift_2, "hh")

def hh_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            distance_label.append(distances[dlist.index(dist_line)])
            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "hh_dist", distance_label)

def chh_dist_3d(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    shift_3 = []  
    shift_4 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for c_ind_1 in c_atm_ind:
                if (simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][c_ind_1] != 0): 
                    for c_ind_2 in c_atm_ind:
                        if (simple_atom_positions[c_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_1][c_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][c_ind_1], atoms[c_ind_1]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
                            shift_3.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][c_ind_2], atoms[c_ind_2]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
                            shift_3.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_3dd(shift_1, shift_2, shift_3, "chh_dist_3d", distance_label)

def nhh_dist_3d(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    shift_3 = []  
    shift_4 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for n_ind_1 in n_atm_ind:
                if (simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][n_ind_1] != 0): 
                    for n_ind_2 in n_atm_ind:
                        if (simple_atom_positions[n_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_1][n_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][n_ind_1], atoms[n_ind_1]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
                            shift_3.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][n_ind_2], atoms[n_ind_2]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
                            shift_3.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_3dd(shift_1, shift_2, shift_3, "nhh_dist_3d", distance_label)

def nhh_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for n_ind_1 in n_atm_ind:
                if (simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][n_ind_1] != 0): 
                    for h_ind_2 in h_atm_ind:
                        if (simple_atom_positions[h_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_2][h_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][n_ind_1], atoms[n_ind_1]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][h_ind_2], atoms[h_ind_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][h_ind_2], atoms[h_ind_2]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][n_ind_1], atoms[n_ind_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "nhh_dist_2d", distance_label)

def nhhc_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for n_ind_1 in n_atm_ind:
                if (simple_atom_positions[n_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][n_ind_1] != 0): 
                    for c_ind_2 in c_atm_ind:
                        if (simple_atom_positions[c_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_2][c_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][n_ind_1], atoms[n_ind_1]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][c_ind_2], atoms[c_ind_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][c_ind_2], atoms[c_ind_2]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][n_ind_1], atoms[n_ind_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "nhhc_dist", distance_label)

def chh_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for c_ind_1 in c_atm_ind:
                if (simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][c_ind_1] != 0): 
                    for h_ind_2 in h_atm_ind:
                        if (simple_atom_positions[h_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_2][h_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][c_ind_1], atoms[c_ind_1]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][h_ind_2], atoms[h_ind_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][h_ind_2], atoms[h_ind_2]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][c_ind_1], atoms[c_ind_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "chh_dist_2d", distance_label)

def chhc_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if (atom_type[pos_1] == 'H') and (atom_type[pos_2] == 'H') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0):
            for c_ind_1 in c_atm_ind:
                if (simple_atom_positions[c_ind_1] == simple_atom_positions[pos_1]) and (chem_shifts[resi_1][c_ind_1] != 0): 
                    for c_ind_2 in c_atm_ind:
                        if (simple_atom_positions[c_ind_2] == simple_atom_positions[pos_2]) and (chem_shifts[resi_2][c_ind_2] != 0):
                            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][c_ind_1], atoms[c_ind_1]])
                            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][c_ind_2], atoms[c_ind_2]])
                            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][c_ind_2], atoms[c_ind_2]])
                            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][c_ind_1], atoms[c_ind_1]])
                            distance_label.append(distances[dlist.index(dist_line)])
                            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "chhc_dist", distance_label)

def cc_spindiff_intra(chem_shifts):
    #this function produces peak list for intra-residue CC spin diffusion experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atom in enumerate(atom_type): #loop through all atomtypes
            if atom == 'C': #check if it is a carbon atom
                for c_ind in c_atm_ind: #loop through all the carbon indices
                    if (c_ind != j) and (chem_shifts[i][c_ind] != 0) and (chem_shifts[i][j]) != 0:
                        shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                        shift_2.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
    write_2d(shift_1, shift_2, "cc_spindiff_intra")

def cc_spindiff_inter(chem_shifts):
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence):
        for j, atom in enumerate(atom_type):
            if atom == 'C':
                for c_ind in c_atm_ind:
                    if (c_ind != j) and (chem_shifts[i][c_ind] != 0) and (chem_shifts[i][j] != 0):
                        shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                        shift_2.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
                if i != 0:
                    for c_ind in c_atm_ind: 
                        if (chem_shifts[i-1][c_ind] != 0) and (chem_shifts[i][j] != 0):
                            shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                            shift_2.append([sequence[i-1], i, chem_shifts[i-1][c_ind], atoms[c_ind]])
                if i+1 != len(sequence):
                    for c_ind in c_atm_ind:
                        if (chem_shifts[i+1][c_ind] != 0 and chem_shifts[i][j] != 0):
                            shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                            shift_2.append([sequence[i+1], i+2, chem_shifts[i+1][c_ind], atoms[c_ind]])
    write_2d(shift_1, shift_2, "cc_spindiff_inter")

def cc_spindiff_dist(chem_shifts, dlist, distances):
    shift_1 = []  
    shift_2 = []  
    distance_label = []  
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if ((atom_type[pos_1] == 'C') and (atom_type[pos_1] == 'C') and (chem_shifts[resi_1][pos_1] != 0) and (chem_shifts[resi_2][pos_2] != 0)):
            shift_1.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            shift_2.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            shift_1.append([sequence[resi_2], (resi_2)+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            shift_2.append([sequence[resi_1], (resi_1)+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            distance_label.append(distances[dlist.index(dist_line)])
            distance_label.append(distances[dlist.index(dist_line)])
    write_2dd(shift_1, shift_2, "cc_spindiff_dist", distance_label)


def nca(chem_shifts):
    #this function produces peak list for intra-residue NCA experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((chem_shifts[i][atoms.index('N')] != 0) 
            and (chem_shifts[i][atoms.index('CA')] != 0)): #check if the chemical shifts of the two are not zero
            shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
            shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
    write_2d(shift_1, shift_2, "nca")


def nco(chem_shifts):
    #this function produces peak list for inter-residue NCO experiment
    shift_1 = []
    shift_2 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((i != 0) 
            and (chem_shifts[i][atoms.index('N')] != 0)
            and (chem_shifts[i-1][atoms.index('C')] != 0)): #check if the chemical shifts of the two are not zero
            shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
            shift_2.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('C')], "C"])
    write_2d(shift_1, shift_2, "nco")

def ncocx_dist(chem_shifts, dimensionality, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if ((atm_1 == 'C') 
            and (chem_shifts[resi_1+1][atoms.index('N')] != 0) 
            and (chem_shifts[resi_1][pos_1] != 0) 
            and (chem_shifts[resi_2][pos_2] != 0)):
            shift_1.append([sequence[(resi_1)+1], (resi_1)+2, chem_shifts[(resi_1)+1][atoms.index('N')], "N"])
            shift_2.append([sequence[resi_1], resi_1+1, chem_shifts[resi_1][pos_1], "CA"])
            shift_3.append([sequence[resi_2], resi_2+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            distance_label.append(distances[dlist.index(dist_line)])
        elif ((atm_2 == 'C') 
              and (chem_shifts[resi_2+1][atoms.index('N')] != 0) 
              and (chem_shifts[resi_2][pos_2] != 0) 
              and (chem_shifts[resi_1][pos_1] != 0)):
            shift_1.append([sequence[(resi_2)+1], (resi_2)+2, chem_shifts[(resi_2)+1][atoms.index('N')], "N"])
            shift_2.append([sequence[resi_2], resi_2+1, chem_shifts[resi_2][pos_2], "CA"])
            shift_3.append([sequence[resi_1], resi_1+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            distance_label.append(distances[dlist.index(dist_line)])
    if dimensionality == 2: write_2dd(shift_1, shift_3, "ncocx_dist_2d", distance_label)
    elif dimensionality == 3: write_3dd(shift_1, shift_2, shift_3, "ncocx_dist_3d", distance_label)

def ncacx_dist(chem_shifts, dimensionality, dlist, distances):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    distance_label = []
    for dist_line in dlist:
        resi_1 = int(dist_line[0]) -1
        atm_1 = dist_line[1]
        resi_2 = int(dist_line[2]) - 1
        atm_2 = dist_line[3]
        pos_1 = atoms.index(atm_1)
        pos_2 = atoms.index(atm_2)
        if ((atm_1 == 'CA') 
            and (chem_shifts[resi_1][atoms.index('N')] != 0) 
            and (chem_shifts[resi_1][pos_1] != 0) 
            and (chem_shifts[resi_2][pos_2] != 0)):
            shift_1.append([sequence[resi_1], resi_1+1, chem_shifts[resi_1][atoms.index('N')], "N"])
            shift_2.append([sequence[resi_1], resi_1+1, chem_shifts[resi_1][pos_1], "CA"])
            shift_3.append([sequence[resi_2], resi_2+1, chem_shifts[resi_2][pos_2], atoms[pos_2]])
            distance_label.append(distances[dlist.index(dist_line)])
        elif ((atm_2 == 'CA') 
              and (chem_shifts[resi_2][atoms.index('N')] != 0) 
              and (chem_shifts[resi_2][pos_2] != 0)
              and (chem_shifts[resi_1][pos_1] != 0)):
            shift_1.append([sequence[resi_2], resi_2+1, chem_shifts[resi_2][atoms.index('N')], "N"])
            shift_2.append([sequence[resi_2], resi_2+1, chem_shifts[resi_2][pos_2], "CA"])
            shift_3.append([sequence[resi_1], resi_1+1, chem_shifts[resi_1][pos_1], atoms[pos_1]])
            distance_label.append(distances[dlist.index(dist_line)])
    if dimensionality == 2: write_2dd(shift_1, shift_3, "ncacx_dist_2d", distance_label)
    elif dimensionality == 3: write_3dd(shift_1, shift_2, shift_3, "ncacx_dist_3d", distance_label)
    
def ncacx_inter(chem_shifts, dimensionality):
    #this function produces peak list for inter residue NCACX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((chem_shifts[i][atoms.index('N')] != 0) 
            and (chem_shifts[i][atoms.index('CA')] != 0)): #check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if ((chem_shifts[i][c_ind] != 0) and (atoms[c_ind] != 'CA')):
                    shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                    shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
                    shift_3.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
            if i != 0:
                for c_ind in c_atm_ind:
                    if ((chem_shifts[i-1][c_ind] != 0) and (atoms[c_ind] != 'CA')):
                        shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                        shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
                        shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind], atoms[c_ind]])
            if i+1 != len(sequence):
                for c_ind in c_atm_ind:
                    if ((chem_shifts[i+1][c_ind] != 0) and (atoms[c_ind] != 'CA')):
                        shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                        shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
                        shift_3.append([sequence[i+1], i+2, chem_shifts[i+1][c_ind], atoms[c_ind]])
    if dimensionality == 2: write_2d(shift_1, shift_3, "ncacx_inter_2d")
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, "ncacx_inter_3d")

def ncacx(chem_shifts, dimensionality):
    #this function produces peak list for inter residue NCACX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((chem_shifts[i][atoms.index('H')] != 0) 
            and (chem_shifts[i][atoms.index('N')] != 0) 
            and (chem_shifts[i][atoms.index('CA')] != 0)): #check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if ((chem_shifts[i][c_ind] != 0) and (atoms[c_ind] != 'CA')):
                    shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                    shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
                    shift_3.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
    if dimensionality == 2: write_2d(shift_1, shift_3, "ncacx_2d")
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, "ncacx_3d")

def ncocx(chem_shifts, dimensionality):
    #this function produces peak list for inter residue NCOCX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((i != 0)
            and (chem_shifts[i][atoms.index('N')] != 0)
            and (chem_shifts[i-1][atoms.index('C')] != 0)): #check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if ((chem_shifts[i-1][c_ind] != 0) and (atoms[c_ind] != 'C')):
                    shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                    shift_2.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('C')], "C"])
                    shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind], atoms[c_ind]])
    if dimensionality == 2: write_2d(shift_1, shift_3, "ncocx_2d")
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, "ncocx_3d")

def ncoca_cb(chem_shifts, dimensionality):
    #this function produces peak list for inter residue NCOCX experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((i != 0)
            and (chem_shifts[i][atoms.index('N')] != 0)
            and (chem_shifts[i-1][atoms.index('C')] != 0)): #check if the chemical shifts of the two are not zero
            for c_ind in c_atm_ind:
                if ((chem_shifts[i-1][c_ind] != 0) and (atoms[c_ind] == 'CA' or atoms[c_ind] == 'CB')):
                    shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                    shift_2.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('C')], "C"])
                    shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind], atoms[c_ind]])
    if dimensionality == 2: write_2d(shift_1, shift_3, "ncocacb_2d")
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, "ncocacb_3d")

def dqsqsq_inter(chem_shifts):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    neighbors = [['C', 'A'], ['A', 'B'], ['B', 'G'], ['G', 'D'], ['D', 'E'], ['E', 'Z'], ['H', 'Z'], ['E', 'H']]
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (atom_type[j] == 'C') and (j+1+k < len(atom_positions)) and ([atm_pos, simple_atom_positions[j+1+k]] in neighbors):
                    if (chem_shifts[i][j] != 0) and (chem_shifts[i][j+1+k] != 0):
                        shift_1.append([residue, i+1, chem_shifts[i][j] + chem_shifts[i][j+1+k], '%s+%s' %(atoms[j], atoms[j+1+k])])
                        shift_2.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                        shift_3.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                        shift_1.append([residue, i+1, chem_shifts[i][j+1+k] + chem_shifts[i][j], '%s+%s' %(atoms[j+1+k], atoms[j])])
                        shift_2.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                        shift_3.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                if i != 0:
                    if (atom_type[j] == 'C') and (j+1+k < len(atom_positions)) and ([atm_pos, simple_atom_positions[j+1+k]] in neighbors):
                        if (chem_shifts[i][j] != 0) and (chem_shifts[i][j+1+k] != 0):
                            for c_ind in c_atm_ind:
                                shift_1.append([residue, i+1, chem_shifts[i][j] + chem_shifts[i][j+1+k], '%s+%s' %(atoms[j], atoms[j+1+k])])
                                shift_2.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                                shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind], '%s' %(atoms[c_ind])])
                                shift_1.append([residue, i+1, chem_shifts[i][j+1+k] + chem_shifts[i][j], '%s+%s' %(atoms[j+1+k], atoms[j])])
                                shift_2.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                                shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind], '%s' %(atoms[c_ind])])
                if i != len(sequence)-1:
                    if (atom_type[j] == 'C') and (j+1+k < len(atom_positions)) and ([atm_pos, simple_atom_positions[j+1+k]] in neighbors):
                        if (chem_shifts[i][j] != 0) and (chem_shifts[i][j+1+k] != 0):
                            for c_ind in c_atm_ind:
                                shift_1.append([residue, i+1, chem_shifts[i][j] + chem_shifts[i][j+1+k], '%s+%s' %(atoms[j], atoms[j+1+k])])
                                shift_2.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                                shift_3.append([sequence[i+1], i+2, chem_shifts[i-1][c_ind], '%s' %(atoms[c_ind])])
                                shift_1.append([residue, i+1, chem_shifts[i][j+1+k] + chem_shifts[i][j], '%s+%s' %(atoms[j+1+k], atoms[j])])
                                shift_2.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                                shift_3.append([sequence[i+1], i+2, chem_shifts[i-1][c_ind], '%s' %(atoms[c_ind])])
    write_2d(shift_1, shift_2, 'dqsqsq_inter')
    
def dqsqsq_intra(chem_shifts):
    shift_1 = []
    shift_2 = []
    shift_3 = []
    neighbors = [['C', 'A'], ['A', 'B'], ['B', 'G'], ['G', 'D'], ['D', 'E'], ['E', 'Z'], ['H', 'Z'], ['E', 'H']]
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (atom_type[j] == 'C') and (j+1+k < len(atom_positions)) and ([atm_pos, simple_atom_positions[j+1+k]] in neighbors):
                    if (chem_shifts[i][j] != 0) and (chem_shifts[i][j+1+k] != 0):
                        shift_1.append([residue, i+1, chem_shifts[i][j] + chem_shifts[i][j+1+k], '%s+%s' %(atoms[j], atoms[j+1+k])])
                        shift_2.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                        shift_3.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                        shift_1.append([residue, i+1, chem_shifts[i][j+1+k] + chem_shifts[i][j], '%s+%s' %(atoms[j+1+k], atoms[j])])
                        shift_2.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
                        shift_3.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
    write_3d(shift_1, shift_2, shift_3, 'dqsqsq_intra')

def dqsq(chem_shifts):
    shift_1 = []
    shift_2 = []
    neighbors = [['C', 'A'], ['A', 'B'], ['B', 'G'], ['G', 'D'], ['D', 'E'], ['E', 'Z'], ['H', 'Z'], ['E', 'H']]
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atm_pos in enumerate(simple_atom_positions):
            for k in range(6):
                if (atom_type[j] == 'C') and (j+1+k < len(atom_positions)) and ([atm_pos, simple_atom_positions[j+1+k]] in neighbors):
                    if (chem_shifts[i][j] != 0) and (chem_shifts[i][j+1+k] != 0):
                        shift_1.append([residue, i+1, chem_shifts[i][j] + chem_shifts[i][j+1+k], '%s+%s' %(atoms[j], atoms[j+1+k])])
                        shift_2.append([residue, i+1, chem_shifts[i][j], '%s' %(atoms[j])])
                        shift_1.append([residue, i+1, chem_shifts[i][j+1+k] + chem_shifts[i][j], '%s+%s' %(atoms[j+1+k], atoms[j])])
                        shift_2.append([residue, i+1, chem_shifts[i][j+1+k], '%s' %(atoms[j+1+k])])
    write_2d(shift_1, shift_2, 'dqsq')

def canh(chem_shifts, dimensionality):
    #this function produces the peak list for proton detected CA-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if ((chem_shifts[i][atoms.index('H')] != 0) 
            and (chem_shifts[i][atoms.index('N')] != 0)
            and (chem_shifts[i][atoms.index('CA')] != 0)): #check if the chemical shifts of the CA, HA and N are not zero
            shift_1.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
            shift_2.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
            shift_3.append([residue, i+1, chem_shifts[i][atoms.index('H')], "H"])
    if dimensionality == 2: write_2d(shift_1, shift_3, 'canh_2d')
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, 'canh_3d')

def conh(chem_shifts, dimensionality):
    #this function produces the peak list for proton detected C-N(n+1)-HA(n+1) experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if i != 0:
            if ((chem_shifts[i-1][atoms.index('H')] != 0)
                and (chem_shifts[i-1][atoms.index('N')] != 0)
                and (chem_shifts[i][atoms.index('C')] != 0)): #check if the chemical shifts of the Co, H and N are not zero
                shift_1.append([residue, i+1, chem_shifts[i][atoms.index('C')], "C"])
                shift_2.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('N')], "N"])
                shift_3.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('H')], "H"])
    if dimensionality == 2: write_2d(shift_1, shift_3, 'conh_2d')
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, 'conh_3d')
        
def caconh(chem_shifts, dimensionality):
    #this function produces the peak list for proton detected CA-(CO)-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if i != 0:
            if ((chem_shifts[i][atoms.index('H')] != 0) 
                and (chem_shifts[i][atoms.index('N')] != 0)
                and (chem_shifts[i-1][atoms.index('CA')] != 0)
                and (chem_shifts[i-1][atoms.index('C')] != 0)): 
                shift_1.append([sequence[i-1], i, chem_shifts[i-1][atoms.index('CA')], "CA"])
                shift_2.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                shift_3.append([residue, i+1, chem_shifts[i][atoms.index('H')], "H"])
    if dimensionality == 2: write_2d(shift_1, shift_3, 'caconh')
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, 'caconh')

def cocanh(chem_shifts, dimensionality):
    #this function produces the peak list for proton detected CO-(CA)-N-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if i != 0:
            if ((chem_shifts[i-1][atoms.index('H')] != 0)
                and (chem_shifts[i-1][atoms.index('N')] != 0)
                and (chem_shifts[i][atoms.index('C')] != 0)
                and (chem_shifts[i][atoms.index('CA')] != 0)): 
                shift_1.append([residue, i+1, chem_shifts[i][42], "C"])
                shift_2.append([sequence[i-1], i, chem_shifts[i-1][33], "N"])
                shift_3.append([sequence[i-1], i, chem_shifts[i-1][0], "H"])
    if dimensionality == 2: write_2d(shift_1, shift_3, 'cocanh_2d')
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, 'cocanh_3d')

def ncah(chem_shifts, dimensionality):
    #this function produces the peak list for proton detected N-CA-H experiment
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        if i != 0:
            if ((chem_shifts[i][atoms.index('H')] != 0)
                and (chem_shifts[i][atoms.index('N')] != 0)
                and (chem_shifts[i][atoms.index('HA')] != 0)): 
                shift_1.append([residue, i+1, chem_shifts[i][atoms.index('N')], "N"])
                shift_2.append([residue, i+1, chem_shifts[i][atoms.index('CA')], "CA"])
                shift_3.append([residue, i+1, chem_shifts[i][atoms.index('HA')], "HA"])
    if dimensionality == 2: write_2d(shift_1, shift_3, 'ncah_2d')
    elif dimensionality == 3: write_3d(shift_1, shift_2, shift_3, 'ncah_3d')

def sqsqsq_inter(chem_shifts):
    #this function produces peak list for intra-residue CC spin diffusion experiment where
    shift_1 = []
    shift_2 = []
    shift_3 = []
    for i, residue in enumerate(sequence): #loop through the sequence
        for j, atom in enumerate(atom_type): #loop through all atomtypes
            if atom == 'C': #check if it is a carbon atom
                for c_ind in c_atm_ind: #loop through all the carbon indices
                    if c_ind != j: #ensure that the two indices "j" and "c_ind" are not the same
                        for c_ind_1 in c_atm_ind: #loop through all the carbon indices for i-1
                            if i != 0:
                                if chem_shifts[i-1][c_ind_1] != 0 and chem_shifts[i][c_ind] != 0 and chem_shifts[i][j] != 0: #check if the chemical shifts of the three are not zero
                                    shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                                    shift_2.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
                                    shift_3.append([sequence[i-1], i, chem_shifts[i-1][c_ind_1], atoms[c_ind_1]])
                        for c_ind_1 in c_atm_ind: #loop through all the carbon indices for i+1
                            if i != len(sequence) -1:
                                if chem_shifts[i+1][c_ind_1] != 0 and chem_shifts[i][c_ind] != 0 and chem_shifts[i][j] != 0: #check if the CS are not zero
                                    shift_1.append([residue, i+1, chem_shifts[i][j], atoms[j]])
                                    shift_2.append([residue, i+1, chem_shifts[i][c_ind], atoms[c_ind]])
                                    shift_3.append([sequence[i+1], i+2, chem_shifts[i+1][c_ind_1], atoms[c_ind_1]])
    write_3d(shift_1, shift_2, shift_3, 'sqsqsq_inter')

def write_2d(shift_a, shift_b, extension):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = "%s%s%s-%s%s%s" %(shift_a[i][0], shift_a[i][1], shift_a[i][3], shift_b[i][0], shift_b[i][1], shift_b[i][3])
        if sl == True: column_1 = column_4
        else: column_1 = "?-?"
        peak.append([column_1, column_2, column_3, column_4])
    with open('%s/%s-%s.txt' %(w_dir, project_name, extension), 'w') as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\n" %(pk[0], pk[1], pk[2], pk[3]))

def write_2dd(shift_a, shift_b, extension, distance_label):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = "%s%s%s-%s%s%s_%s" %(shift_a[i][0], shift_a[i][1], shift_a[i][3], shift_b[i][0], shift_b[i][1], shift_b[i][3], distance_label[i])
        if sl == True: column_1 = column_4
        else: column_1 = "?-?"
        peak.append([column_1, column_2, column_3, column_4])
    with open('%s/%s-%s.txt' %(w_dir, project_name, extension), 'w') as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\n" %(pk[0], pk[1], pk[2], pk[3]))

def write_3dd(shift_a, shift_b, shift_c, extension, distance_label):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = round(shift_c[i][2], 4)
        column_5 = "%s%s%s-%s%s%s-%s%s%s_%s" %(shift_a[i][0], shift_a[i][1], shift_a[i][3], 
                                               shift_b[i][0], shift_b[i][1], shift_b[i][3],
                                               shift_c[i][0], shift_c[i][1], shift_c[i][3], distance_label[i])
        if sl == True: column_1 = column_5
        else: column_1 = "?-?-?"
        peak.append([column_1, column_2, column_3, column_4, column_5])
    with open('%s/%s-%s.txt' %(w_dir, project_name, extension.lower()), 'w') as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\t%s\n" %(pk[0], pk[1], pk[2], pk[3], pk[4]))

def write_3d(shift_a, shift_b, shift_c, extension):
    peak = []
    for i in range(len(shift_a)):
        column_2 = round(shift_a[i][2], 4)
        column_3 = round(shift_b[i][2], 4)
        column_4 = round(shift_c[i][2], 4)
        column_5 = "%s%s%s-%s%s%s-%s%s%s" %(shift_a[i][0], shift_a[i][1], shift_a[i][3], 
                                            shift_b[i][0], shift_b[i][1], shift_b[i][3],
                                            shift_c[i][0], shift_c[i][1], shift_c[i][3])
        if sl == True: column_1 = column_5
        else: column_1 = "?-?-?"
        peak.append([column_1, column_2, column_3, column_4, column_5])
    with open('%s/%s-%s.txt' %(w_dir, project_name, extension.lower()), 'w') as output:
        for pk in peak:
            output.write("%s\t%s\t%s\t%s\t%s\n" %(pk[0], pk[1], pk[2], pk[3], pk[4]))
    
def fractional_deuteration(sequence, chem_shifts):
    for i, residue in enumerate(sequence):
        if residue == 'A': rem_atm = ['HA']
        elif residue == 'R': rem_atm = ['HA', 'HB']
        elif residue == 'D': rem_atm = ['HA']
        elif residue == 'N': rem_atm = ['HA']
        elif residue == 'C': rem_atm = ['HA']
        elif residue == 'E': rem_atm = ['HA', 'HB']
        elif residue == 'Q': rem_atm = ['HA', 'HB']
        elif residue == 'G': rem_atm = ['HA']
        elif residue == 'H': rem_atm = ['HA']
        elif residue == 'I': rem_atm = ['HA', 'HB', 'HG']
        elif residue == 'K': rem_atm = ['HA']
        elif residue == 'M': rem_atm = ['HA']
        elif residue == 'P': rem_atm = ['HA', 'HB']
        elif residue == 'L': rem_atm = ['HA', 'HB', 'HG']
        elif residue == 'F': rem_atm = ['HA']
        elif residue == 'S': rem_atm = ['HA']
        elif residue == 'T': rem_atm = ['HA']
        elif residue == 'W': rem_atm = ['HA']
        elif residue == 'Y': rem_atm = ['HA']
        elif residue == 'V': rem_atm = ['HA', 'HB']
        for j, atom in enumerate(atoms):
            for element in rem_atm:
                if atom == element: chem_shifts[i][j] = 0.00
    return(chem_shifts)
    
def glycerol_label(sequence, chem_shifts, label_type):
    if label_type == 13:
        for i, residue in enumerate(sequence):
            if residue == 'A': rem_atm = ['CA', 'CG1', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'R': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'D': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'N': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'C': rem_atm = ['CA', 'CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'E': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'Q': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'G': rem_atm = ['CA', 'CB','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'H': rem_atm = ['CA', 'CG2', 'CD1', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'I': rem_atm = ['CB', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'L': rem_atm = ['C', 'CB','CG1','CG2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'K': rem_atm = ['CG2', 'CD2', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'M': rem_atm = ['CG2', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'F': rem_atm = ['CA', 'CG1','CG2', 'CE2', 'CE3', 'CZ2', 'CZ3','CH']
            elif residue == 'P': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'S': rem_atm = ['CA','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'T': rem_atm = ['CG1', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'W': rem_atm = ['CA', 'CG2', 'CD2', 'CE1', 'CZ1', 'CZ3']
            elif residue == 'Y': rem_atm = ['CA', 'CG1','CG2', 'CE2', 'CE3', 'CZ2', 'CZ3','CH']
            elif residue == 'V': rem_atm = ['CA', 'CB', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            for j, atom in enumerate(atoms):
                for element in rem_atm:
                    if atom == element: chem_shifts[i][j] = 0.00
    elif label_type == 2:
        for i, residue in enumerate(sequence):
            if residue == 'A': rem_atm = ['C', 'CB','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'R': rem_atm = ['CG1','CG2', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'D': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'N': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'C': rem_atm = ['C', 'CB','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'E': rem_atm = ['CG1','CG2', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'Q': rem_atm = ['CG1','CG2', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'G': rem_atm = ['C', 'CB','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'H': rem_atm = ['C', 'CB', 'CG2', 'CD1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'I': rem_atm = ['CG2', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'L': rem_atm = ['CA', 'CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'K': rem_atm = ['CG2', 'CD2', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'M': rem_atm = ['CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'F': rem_atm = ['C', 'CB', 'CG2', 'CD1', 'CD2', 'CE1', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'P': rem_atm = ['CG1','CG2', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'S': rem_atm = ['C', 'CB','CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'T': rem_atm = ['CG1', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'W': rem_atm = ['C', 'CB', 'CG2', 'CD1', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CH']
            elif residue == 'Y': rem_atm = ['C', 'CB', 'CG2', 'CD1', 'CD2', 'CE1', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            elif residue == 'V': rem_atm = ['C', 'CG1','CG2', 'CD1', 'CD2', 'CE1', 'CE2', 'CE3', 'CZ1', 'CZ2', 'CZ3','CH']
            for j, atom in enumerate(atoms):
                for element in rem_atm:
                    if atom == element: chem_shifts[i][j] = 0.00
    else:
        print "ERROR: %s is an invalid glycerol labelling scheme" %label_type
        sys.exit()
    return(chem_shifts)

def rev_label(chem_shifts, order, amino_acid):
    if order == 'rev_cn':
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for n_ind in n_atm_ind:
                    chem_shift[n_ind] = 0.00
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for c_ind in c_atm_ind:
                    chem_shift[c_ind] = 0.00
    if order == 'rev_n':
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for n_ind in n_atm_ind:
                    chem_shift[n_ind] = 0.00
    if order == 'rev_c':
        for i, chem_shift in enumerate(chem_shifts):
            if sequence[i] == amino_acid:
                for c_ind in c_atm_ind:
                    chem_shift[c_ind] = 0.00
    return(chem_shifts)
        
def replace_bmrb(chem_shifts, bmrb_tables_file, bmrb_columns):
    res_num = int(bmrb_columns[0]) - 1 #Column index for residue number
    atm_nam =  int(bmrb_columns[1]) - 1 #Column index for atom name
    c_s =  int(bmrb_columns[2]) - 1 #Column index for chemical shift
    bmrb_table = []
    with open("%s/%s" %(w_dir, bmrb_tables_file), 'r') as bmrb_file:
        for line in bmrb_file:
            bmrb_table.append(line.split())
        bmrb_file.close()
    for c_shift in bmrb_table:
        for i, atom in enumerate(atoms):
            if c_shift[atm_nam] == atom:
                 try: 
                     res_ind = int(c_shift[res_num]) - 1
                     chem_shifts[res_ind][i] = float(c_shift[c_s])
                 except:
                     continue
    return(chem_shifts)

def assign_chemical_shifts(sequence, sec_struc):
    reference = np.fromfile("%s/standard.dat" %w_dir, dtype = float, sep = " ").reshape(((4, 20, 59)))
    shifts = np.zeros((len(sequence), 59))
    for i, residue in enumerate(sequence):
        for j, amino_acid in enumerate(amino_acids):
            if residue == amino_acid:
                for k, secondary_structure in enumerate(secondary_structures):
                    if sec_struc[i] == secondary_structure:
                        shifts[i] = reference[k, j]
    return(shifts)

def check_user_input(user_input, input_type, error_type):
    matches = 0
    for component in user_input:
        for reference_component in input_type:
            if component == reference_component: 
                matches += 1
    if matches != len(user_input):
        print "%s Error: Please check!" %error_type
        sys.exit()

if __name__ == '__main__':
    main()
