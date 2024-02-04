import os
import sys

# Open the file containing PDB IDs
data = open('test.txt', 'r')
ids = data.readlines()

# Dictionary mapping three-letter amino acid codes to one-letter codes
aa_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
           'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
           'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

# Specify the directory for the dataset
datasetDir = 'pdbDataset_04022024'

# Check if the directory already exists, and create it if not
if not os.path.exists(datasetDir + '/SHEET'):
    os.mkdir(datasetDir + '/SHEET')

# Iterate through each PDB ID
for ID in ids:
    n = ID[0:4]

    # Open the PDB file for reading
    in_file = open(datasetDir + '/' + n + '.pdb', 'r')
    
    # Open the output file for writing SHEET information
    out_file = open(datasetDir + '/SHEET/' + n + '.pdb', 'a+')
    out_file.write('START\n')

    # Read all lines from the input file
    infile_lines = in_file.readlines()
    helix_initial = list()
    sheet_initial = list()
    sheet_terminal = list()
    helix_terminal = list()

    # Iterate through each line in the PDB file
    for line in infile_lines:
        if line.startswith('HELIX'):
            if line[19:20].strip() == 'A':
                ini = int(line[22:25].strip())
                ter = int(line[34:37].strip())
                helix_initial.append(ini)
                helix_terminal.append(ter)

        if line.startswith('SHEET'):
            if line[21:22].strip() == 'A':
                ini = int(line[23:26].strip())
                ter = int(line[34:37].strip())
                sheet_initial.append(ini)
                sheet_terminal.append(ter)

    # Extract SHEET information
    for l in range(len(sheet_initial)):
        end = sheet_terminal[l]
        seq = ""
        for i in range(sheet_initial[l] - 1, end + 2):
            for line in infile_lines:
                if line.startswith('ATOM'):
                    if line[21:22].strip() == 'A':
                        if int(line[23:26].strip()) == int(i):
                            out_file.write(line)
                            if line[13:15] == 'CA':
                                res = str(line[17:20])
                                seq = seq + aa_dict[res]
        out_file.write("SEQ\t" + str(seq) + "\nEND\nSTART\n")
        del seq
        
'''
        #for extracting just HELIX 
        
        for l in range(len(helix_initial)):
                end=helix_terminal[l]
                seq=""
                for i in range(helix_initial[l],end+1):
                        
                        for line in infile_lines:
                                if line.startswith('ATOM'):
                                        if line[21:22].strip() == 'A':
                                                if int(line[23:26].strip()) == int(i):
                                                        out_file.write(line)
                                                        if line[13:15] == 'CA':
                                                                res=str(line[17:20])
                                                                seq=seq+aa_dict[res]
                out_file.write("SEQ\t"+str(seq)+"\nEND\nSTART\n")
                del seq

             
        #for extracting loops between two helices
       
        hlh_initial=list()
        hlh_terminal=list()
        
        for i in range(0,len(helix_initial)):
                try:
                        count=0
                        h_t=int(helix_terminal[i])
                
                        h_i=int(helix_initial[i+1])
                        for s in range(len(sheet_initial)):
                                s_i=sheet_initial[s]
                                if int(h_i) > int(s_i) > int(h_t):
                                        count=count+1
                        if count == 0:
                               # one_letter=list()
                                seq= ''
                                for i in range(int(h_t-4),int(h_i+5)):
                                        for line in infile_lines:
                                                if line.startswith('ATOM'):
                                                        if line[21:22].strip() == 'A':
                                                                if int(line[23:26].strip()) == int(i):
                                                                        out_file.write(line)
                                           
                                                                        if line[13:15] == 'CA':
                                                                                res=str(line[17:20])
                                                                                seq=seq+aa_dict[res]
                                out_file.write("SEQ\t"+str(seq)+"\nEND\nSTART\n")
                                del seq
                                
                except:
                        break
       
       #for extracting loops between two SHEETS
     
        sheet_initial=sorted(sheet_initial)
        sheet_terminal=sorted(sheet_terminal)
        for i in range(0,len(sheet_initial)-1):
                try:
                        count=0
                        s_t=int(sheet_terminal[i])
                
                        s_i=int(sheet_initial[i+1])
                        for s in range(len(helix_initial)):
                                h_i=helix_initial[s]
                                if int(s_i) > int(h_i) > int(s_t):
                                                                       
                                        count=count+1
                        if count == 0:
                                seq=""
                                for i in range(int(s_t-4),int(s_i+5)):
                                        for line in infile_lines:
                                                if line.startswith('ATOM'):
                                                        if line[21:22].strip() == 'A':
                                                                if int(line[23:26].strip()) == int(i):
                                                                        out_file.write(line)
                                                                                
                                                                        if line[13:15] == 'CA':

                                                                                res=str(line[17:20])
                                                                                seq=seq+aa_dict[res]
                                out_file.write("SEQ\t"+str(seq)+"\nEND\nSTART\n")
                                del seq
                except:
                        break
        #for extracting loops between HELIX and SHEET
     
        sheet_initial=sorted(sheet_initial)
        for i in range(0,len(helix_initial)):
                try:
                        count=0
                        h_t=int(helix_terminal[i])
                
                        h_i=int(helix_initial[i+1])
                        for s in range(len(sheet_initial)):
                                s_i=sheet_initial[s]
                                if int(h_i) > int(s_i) > int(h_t):
                                                                       
                                        count=count+1
                                        break
                        if count > 0:
                                seq=''
                                for i in range(int(h_t-4),int(s_i)+5):
                                        for line in infile_lines:
                                                if line.startswith('ATOM'):
                                                        if line[21:22].strip() == 'A':
                                                                if int(line[23:26].strip()) == int(i):
                                                                        out_file.write(line)
                                                                        if line[13:15] == 'CA':

                                                                                res=str(line[17:20])
                                                                                seq=seq+aa_dict[res]
                                out_file.write("SEQ\t"+str(seq)+"\nEND\nSTART\n")
                                del seq
                except:
                        break

        #for extracting loops between SHEET and HELIX
        
     
        sheet_initial=sorted(sheet_initial)
        sheet_terminal=sorted(sheet_terminal)
        for i in range(0,len(sheet_initial)):
                try:
                        count=0
                        s_t=int(sheet_terminal[i])
                
                        s_i=int(sheet_initial[i+1])
                        for h in range(len(helix_initial)):
                                h_i=helix_initial[h]
                                if int(s_i) > int(h_i) > int(s_t):
                                                                       
                                        count=count+1
                                        break
                        if count > 0:
                                seq=''
                                for i in range(int(s_t-4),int(h_i+5)):
                                        for line in infile_lines:
                                                if line.startswith('ATOM'):
                                                        if line[21:22].strip() == 'A':
                                                                if int(line[23:26].strip()) == int(i):
                                                                        out_file.write(line)
                                                                        if line[13:15] == 'CA':
                                                                                res=str(line[17:20])
                                                                                seq=seq+aa_dict[res]
                                out_file.write("SEQ\t"+str(seq)+"\nEND\nSTART\n")
                                del seq
                except:
                        break
'''
in_file.close()
out_file.close()





        
                                      
