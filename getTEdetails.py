'''
Code prepares the following output:
TE details for each TE in the genome
1. percentage of types of TEs
2. DNA
3. LTR
4. LINE
5. SINE
6. Unclassified
7. Total TE (bases masked)

Code takes in the .tbl file output from RepeatModeler.
'''

import sys
import os

# Input: .tbl file from RepeatModeler
# Output: .csv file with TE details
def getTEdetails(inputFile, outputFile, fish_dict = None, header = True):
    ''' Function gets TE details from .tbl file and writes to .csv file
    Input: .tbl file from RepeatModeler
    Output: .csv file with TE details
    '''
    # -- Open input filehttps://prod.liveshare.vsengsaas.visualstudio.com/join?60934BB25B2134CC679B727A685407570E13
    out_dict = {}
    out_str = ""
    # the input file is divided into two parts, a preamble with file details
    # such as name, total TE count, the coverage of the genome, etc.
    # and the actual TE details
    # the second part of the file has a human-readible format of a table
    # from this second part, we need the following details: 
    # DNA, LTR, LINE, SINE, and Unclassified
    # print ("Processing file: ", inputFile)
    with open(inputFile, 'r') as inFile:
        # Open output file
        preamble_flag = True # flag to indicate if we are in the preamble
        preamble_marker = 0
        for line in inFile:
            if line.startswith("===================="):
                # first line of the preamble or last line of preamble
                preamble_marker += 1
                if preamble_marker == 2:
                    # last line of preamble
                    preamble_flag = False
                    continue
            if preamble_flag: # total TE
                line = line.strip()
                if line.startswith("file name: "): # file name
                    val = line.split()[2].split(".")[0]
                    out_dict["file"] = val
                # if line.startswith("total"): # total length
                #     val = int(line.split()[2])
                #     out_dict["total length"] = val
                #     continue
                # elif line.startswith("GC"): # GC level
                #     val = float(line.split()[2])
                #     out_dict["GC"] = val
                #     continue
                if line.startswith("bases masked"): # bases masked
                    out_dict["total_TE"] = line.split()[5]
                else:
                    continue
            else: # TE details
                line = line.strip()
                line = line.split()
                if len(line) == 0:
                    continue
                if line[0].startswith("SINEs"):
                    out_dict["SINEs"] = line[4]
                    continue
                elif line[0].startswith("LINEs"):
                    out_dict["LINEs"] = line[4]
                    continue
                elif line[0].startswith("LTR"):
                    out_dict["LTR"] = line[5]
                    continue
                elif line[0].startswith("DNA"):
                    out_dict["DNA"] = line[5]
                    continue
                elif line[0].startswith("Unclassified"):
                    out_dict["Unclassified"] = line[4]
                    continue
                else:
                    continue
    # -- Write to output file
    with open(outputFile, 'a') as outFile:
        if header:
            outFile.write("File,DNA,LTR,LINE,SINE,Unclassified,Total TE\n")
        if (fish_dict is not None) and (out_dict["file"] in fish_dict):
            out_str = fish_dict[out_dict["file"]]["SPECIES"]
        else:
            print( out_dict["file"])
            out_str = out_dict["file"]
        out_str += "," + out_dict["DNA"] + "," + out_dict["LTR"] + "," 
        out_str += out_dict["LINEs"] + "," + out_dict["SINEs"] + "," 
        out_str += out_dict["Unclassified"] + "," + out_dict["total_TE"] + "\n"

        outFile.write(out_str)
    
    # -- return
    return

def get_fishname_from_fishtable(fish_table_filename):
    ''' Function gets fish name from fish table
    '''
    # -- Read the file
    ''' Header of the tab separated file contains the following:
    ALIAS FAMILY ORDER SPECIES COMMON_NAME
    '''
    # In the first line of the file, create the dictionary with all the names
    # from the next line of the line, fill up the dictionary
    # return the dictionary

    first_line_flag = True
    with open(fish_table_filename, 'r') as fi:
        fish_dict = {}
        fish_params = []
        for line in fi:
            line = line.strip()
            line = line.split("\t")
            if first_line_flag:
                # populate the dictionary
                fish_params = line
                first_line_flag = False
                continue
            if len(line)>0:
                fish_dict[line[0]] = {}
                for i in range(1, len(line)):
                    fish_dict[line[0]][fish_params[i]] = line[i]
    return fish_dict


def read_files_dir(directory):
    ''' Function reads all files in a directory and returns a list of file names
    Input: directory
    Output: list of file names
    '''
    files = []
    for filename in os.listdir(directory):
        if filename.endswith(".tbl"):
            files.append(filename)
    return files

def convert_tbl_to_csv(directory, outputfile, fish_dict_file=None):
    ''' Function converts all .tbl files in a directory to .csv files
    Input: directory
    Output: .csv files
    '''
    files = read_files_dir(directory)
    fish_dict = get_fishname_from_fishtable(fish_dict_file) if fish_dict_file is not None else None
    first_tbl = True
    for file in files:
        # -- get full file path
        filepath = os.path.join(directory, file)
        print ("Processing file: ", filepath)
        if first_tbl:
            getTEdetails(filepath, outputfile, fish_dict=fish_dict, header = True)
            first_tbl = False
        getTEdetails(filepath, outputfile, fish_dict=fish_dict, header = False)
    return

if __name__ == '__main__':
    # -- Get input file
    input_dir = sys.argv[1]
    # -- Get output file
    output_file = sys.argv[2]
    # -- Get fish dictionary file
    fish_dict_file = sys.argv[3]
    # -- Convert all .tbl files in directory to .csv
    convert_tbl_to_csv(input_dir, output_file, fish_dict_file)

    # Get TE details
    # getTEdetails(inputFile, outputFile)
    # getTEdetails("lepisosteus_osseus_noadapter_nodups_noemptylines.fasta.tbl", "out_TEdetails.csv")
    # Sample:
    # python3 getTEDetails.py lepisosteus_osseus_noadapter_nodups_noemptylines.fasta.tbl out_TEdetails.csv