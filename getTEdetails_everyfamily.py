'''
Code prepares the out file for the getTEdetails_everyfamily.py
Input: .out file from RepeatModeler
Output: .csv file with TE details per TE subtype and family
        header would be repeatclass/family for each family
        The rows would be the sample names
        The numbers are the percentages of the TE in the sample
'''

import sys
import os
import re
from collections import defaultdict
from tqdm import tqdm

def get_genome_size(inputTBLfile:str)->int:
    ''' Function gets the genome size from the .out file
    Input: .out file from RepeatModeler
    Output: genome size
    '''
    # filename = inputTBLfile.split(".")[0]
    genome_size = 0
    with open(inputTBLfile, 'r') as fi:
        for line in fi:
            line = line.strip()
            if line.startswith("total length:"):
                genome_size = int(line.split()[2])
                break
    return genome_size

def getheader(line):
    ''' Function gets the header details from the .out file
    Input: header line from the .out file
    Output: class_index, position_index
    '''
    line = line.strip() # remove trailing spaces
    line = line.split() # convert to list
    line = [x.strip() for x in line] # remove the spaces in each element in the list
    class_index = line.index("class/family") # get the index of the class/family
    position_index = (line.index("begin"), line.index("end")) # get the index of the begin and end
    return class_index, position_index # return the index of the class/family and begin and end

def getTEdetails(inputRMFile, header = True)->dict:
    line_count = 0
    class_index = 10
    position_index = (5,6)
    fish_dict = defaultdict(lambda: 0)
    with open(inputRMFile, 'r') as inFile:
        for line in inFile:
            line_count += 1
            # if line_count == 2: #get header details
            #     class_index, position_index = getheader(line)
            if line_count <= 3: #ignore the header and the space
                continue
            line = line.strip()
            # -- ignore the line if it ends with *
            if line.endswith("*"):
                continue
            # -- parse the line
            line = line.split()
            line = [x.strip() for x in line]
            # -- get the class and family
            if line[class_index].startswith("SINE") or line[class_index].startswith("LINE") or line[class_index].startswith("LTR") or line[class_index].startswith("DNA") and not line[class_index].endswith("?"):
                class_family = line[class_index] # get the class
            else:
                continue
            # -- get the begin and end
            begin = int(line[position_index[0]])
            end = int(line[position_index[1]])
            # -- get the length
            length = end - begin + 1
            # -- add it to the dictionary
            fish_dict[class_family] += length
    # -- return the fish_dict
    fish_dict = dict(fish_dict)
    return fish_dict

def match_tbl_rm_files(tbl_file_list, rm_file_list):
    ''' Matches a .tbl file with a .out file
    '''
    tbl_rm_pairs = []
    for tbl_file in tbl_file_list:
        tbl_fish = tbl_file.split(".")[0]
        for rm_file in rm_file_list:
            rm_fish = rm_file.split(".")[0]
            if tbl_fish == rm_fish:
                tbl_rm_pairs.append((tbl_file, rm_file))
    return tbl_rm_pairs

def get_fish_file_to_name_dict(namesFile):
    ''' Function gets the fish file to fish name dictionary
    Input: names file
    Output: fish file to fish name dictionary
    '''
    # -- create a dictionary from the namesFile
    fish_file_to_name_dict = {}
    with open(namesFile, 'r') as fi:
        for line in fi:
            line = line.strip()
            line = line.split("\t")
            fish_file_to_name_dict[line[0]] = line[3]
    return fish_file_to_name_dict

def calculate_percentage(fish_dict:dict, genome_size:int)->dict:
    ''' Function calculates the percentage of the TE in the genome
    Input: fish_dict, genome_size
    Output: fish_dict with the percentage
    '''
    for key, value in fish_dict.items():
        fish_dict[key] = round((value/genome_size)*100,3)
    return fish_dict

def write_output_file(tbl_rm_pairs:list, inputTBLDir:str, inputRMDir:str, outputFile:str, fish_file_to_name_dict:dict)->None:
    ''' Function creates a table where rows are the samples and columns are the TE classes
    Input: tbl_rm_pairs, inputTBLDir, inputRMDir, outputFile, fish_file_to_name_dict
    Output: outputFile
    '''
    output_table_dict = {}
    TE_names_list = []
    # -- get the file name by appending the inputTBLDir and inputRMDir to the files
    for pair in tqdm(tbl_rm_pairs):
        tbl_file = os.path.join(inputTBLDir, pair[0])
        rm_out_file = os.path.join(inputRMDir, pair[1])
        # -- get the length of the genome
        genome_size = get_genome_size(tbl_file)
        # -- get the TE details
        fish_dict = getTEdetails(rm_out_file)
        # -- get the fish name
        fish_basic_name = pair[0].split(".")[0]
        fish_name = fish_basic_name if fish_basic_name not in fish_file_to_name_dict else fish_file_to_name_dict[fish_basic_name]
        # -- calculate the percentages
        fish_dict = calculate_percentage(fish_dict, genome_size)
        # -- add the keys to the TE_names_list
        for keys in fish_dict.keys():
            TE_names_list.append(keys)
        # -- add to dictionary
        output_table_dict[fish_name] = fish_dict
    
    # -- convert the TE_names_list to a set to get unique values
    TE_names_list = list(set(TE_names_list))
    # -- print the data in the output file
    with open(outputFile, 'w') as fo:
        # -- print the TE_names_list
        fo.write("Fish\t" + "\t".join(TE_names_list) + "\n")
        # -- print the data
        for key, value in output_table_dict.items():
            fo.write(re.sub(" ", "_",key) + "\t")
            for TE in TE_names_list:
                if TE in value:
                    fo.write(str(value[TE]) + "\t")
                else:
                    fo.write("0\t")
            fo.write("\n")
    return

def main(inputTBLDir, inputRMDir, namesFile, outputFile):
    # -- get the list of files in the TBL directory
    tbl_file_list = os.listdir(inputTBLDir)
    # -- get the list of files in the RMOut directory
    rm_file_list = os.listdir(inputRMDir)
    # -- print the file lists (first 10 to check)
    # print (".tbl files")
    # print(tbl_file_list[:10])
    # print (".out files")
    # print(rm_file_list[:10])

    # -- match the file names
    tbl_rm_pairs = match_tbl_rm_files(tbl_file_list, rm_file_list)
    # print ("Matched .tbl and .out files")
    # print (tbl_rm_pairs[:10])

    # -- get the fish file to fish name dictionary
    fish_file_to_name_dict = get_fish_file_to_name_dict(namesFile)
    # print ("Fish file to fish name dictionary")
    # # -- print the first 10 elements
    # print (list(fish_file_to_name_dict.items())[:10])

    # -- write the output file
    write_output_file(tbl_rm_pairs, inputTBLDir, inputRMDir, outputFile, fish_file_to_name_dict)
    return

if __name__=='__main__':
    # -- get the input for directories with .tbl files and .out files
    inputTBLDir = sys.argv[1]
    inputRMDir = sys.argv[2]
    # -- get the names text file which maps 'fish_xyz' to 'fish_name'
    namesFile = sys.argv[3]
    # -- get the output file name to output the TEName details table
    outputFile = sys.argv[4]
    # -- main function
    main(inputTBLDir, inputRMDir, namesFile, outputFile)
    # python3 getTEdetails_everyfamily.py 'RepeatMaskTables/' 'RepeatMaskOutFiles/' 'names.tab' 'TEDetails_everyfamily.txt'


##-- END OF FILE --##
## -- how to run -- ##
# python3 getTEdetails_everyfamily.py 'RepeatMaskTables/' 'RepeatMaskOutFiles/' 'names.tab' 'TEDetails_everyfamily.txt'