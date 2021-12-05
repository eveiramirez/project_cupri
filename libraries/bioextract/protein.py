"""
MODULE NAME
        protein

VERSION
        [0.0.1]

PYTHON VERSION
        3.9

AUTHOR
        Ignacio Emmanuel Ramirez Bernabe
        Diego
        Melissa

CONTACT
        iramirez@lcg.unam.mx

DESCRIPTION
        Modulo para la obtencion de informacion de proteinas.

CATEGORY
        Protein
"""

from bioanalyser import (get_taxid)
mport pandas as pd
import numpy as num
import os
from difflib import SequenceMatcher
import argparse
import multiprocessing as mp
import time

# parsing the arguments and warning message
ap = argparse.ArgumentParser(description="Create csv files for different gene names in each Orthologues file")
ap.add_argument("-path",required=False,
	help="Path of folder that containing all folders that have all Orthologues")
ap.add_argument("-outpath",required=False,
	help="Path of directory you want to put output files")

args = vars(ap.parse_args())
Path1=args["path"]
Path2=args["outpath"]

Path1 = ""
Path2 = ""

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█'):    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = '\r')
    # Print New Line on Complete
    if iteration == total: 
        print()

#make output dir
def makedir (path):
    if not os.path.exists(path + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name"): 
        os.mkdir(path + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name")
    if not os.path.exists(path + "/" + "Proccessed_Data" + "/" + "Orthologues_with_same_name"):
        os.mkdir(path + "/" + "Proccessed_Data" + "/" +"Orthologues_with_same_name")


#strip "-" from lists
def comparable(l):
    for x in range(len(l)):
        l[x] = l[x].strip("-")
    return(l)


def abs_list(path):
    temp_list = os.listdir(path)
    main_list = list()
    for t in temp_list:
        temp_list_1 = path + "/" + t
        main_list.append(temp_list_1)
    return(main_list)


#function to decide whether the gene name can be consider the same
# the decision is based on dissimility
def same_filter(first, second):                    
    #caluate similarity of each name
    def similar(str1, str2):
        return SequenceMatcher(None, str1, str2).ratio()
    animal_code = ["ECA","SSC","BTA","CHI","MT","OAR","GGA"] #animal code 
    fir_name_1 = list(first)
    sec_name_1 = list(second)
    sim_score = similar(first,second)
    if first == second:
        dec = "Yes"
    elif first != second:
        if sim_score < 0.5:
            dec = "No"
        elif sim_score >= 0.5:            
            if ("-" in fir_name_1) or ("-" in sec_name_1):
                first_split = first.split("-")
                second_split = second.split("-")                
                for code in animal_code:                    
                    if code in first_split:
                        first_split.remove(code)                        
                        first_back = "".join(first_split)                        
                        break
                    else:
                        first_back = "".join(first_split) 
                for code in animal_code:
                    if code in second_split:
                        second_split.remove(code)
                        second_back = "".join(second_split)
                        break 
                    else:
                        second_back = "".join(second_split)
                if first_back == second_back:
                    dec = "Yes"
                elif first_back == second:
                    dec = "Yes"
                elif first == second_back:
                    dec = "Yes"
                else:
                    dec = "No"               
            else:
                dec = "No"
    return(dec)


#compare gene names in each gene and save in dataframe           
def same_gene(num):
    folder = os.listdir(Path1)
    i = 0
    for folder_name in folder:
        file = os.listdir(Path1 + "/" + folder_name)
        #print("Processing files under" + " " + folder_name + " " + "folder")
        i +=1
        abs_path = abs_list(Path1 + "/" + folder_name)
        
        #for name in file:
        #df = pd.read_csv(Path1 + "/" + folder_name + "/" + name) 
        name = os.path.basename(abs_path[num])
        df = pd.read_csv(abs_path[num])
        gene_name = df["Gene name"]
        name_1 = name.strip(folder_name)
        name_2 = name_1.strip("_")
        name_3 = name_2.strip("txt")
        name_4 = name_3.strip(".")
        other_name = df[name_4 + " " + "gene name"]
        gene_name_new = gene_name.str.upper()
        comparable(gene_name_new)
        other_name_new = other_name.str.upper()
        comparable(other_name_new)
        #sub_df = pd.DataFrame()
        #same_df = pd.DataFrame()           
        #printProgressBar(i - 1 , len(folder), prefix = 'Progress:', suffix = 'Complete', length = 50)            
        if folder_name == "Human": #specific comparasion for human
            print("Processing {}".format(name))
            diff_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name" + "/" + name, "w+")
            same_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_same_name" + "/" + name, "w+")
            column_names = list(df.columns.values)
            diff_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            same_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            for term in range(len(gene_name_new)):
                if gene_name_new[term] == "NONE":
                     pass 
                elif other_name_new[term] == "NONE":
                    pass   
                else: 
                    decision = same_filter(other_name_new[term],gene_name_new[term])
                    if decision == "No":
                        ID_list = df.loc[[term]].values
                        ID_list = list(ID_list[0])
                        diff_file.write("{},{},{},{}\n".format(ID_list[0],ID_list[1],ID_list[2],ID_list[3]))
                        #temp_df = pd.DataFrame(list(df.loc[[term]].values), list(df.columns.values)) #change column names of dataframe
                        #sub_df = sub_df.append(temp_df, ignore_index=True) #output the different gene names in to csv files
                    else:
                        ID_list_1 = list(df.loc[[term]].values)
                        ID_list_1 = list(ID_list_1[0])
                        same_file.write("{},{},{},{}\n".format(ID_list_1[0],ID_list_1[1],ID_list_1[2],ID_list_1[3]))
                        #temp_df_1 = pd.DataFrame(list(df.loc[[term]].values))
                        #same_df = same_df.append(temp_df_1)
            diff_file.close()
            same_file.close()
            print("{} finished".format(name))                          
        elif folder_name == "Mouse": #specific comparasion for mouse
            print("Processing {}".format(name))
            diff_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name" + "/" + name, "w+")
            same_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_same_name" + "/" + name, "w+")
            column_names = list(df.columns.values)
            diff_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            same_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            for term in range(len(gene_name_new)):
                if gene_name_new[term] == "NONE":
                    pass 
                elif other_name_new[term] == "NONE":
                    pass   
                else: 
                    decision = same_filter(other_name_new[term],gene_name_new[term])
                    if decision == "No":
                        ID_list = df.loc[[term]].values
                        ID_list = list(ID_list[0])
                        diff_file.write("{},{},{},{}\n".format(ID_list[0],ID_list[1],ID_list[2],ID_list[3]))
                        #temp_df = pd.DataFrame(list(df.loc[[term]].values), list(df.columns.values)) #change column names of dataframe
                        #sub_df = sub_df.append(temp_df, ignore_index=True) #output the different gene names in to csv files
                    else:
                        ID_list_1 = list(df.loc[[term]].values)
                        ID_list_1 = list(ID_list_1[0])
                        same_file.write("{},{},{},{}\n".format(ID_list_1[0],ID_list_1[1],ID_list_1[2],ID_list_1[3]))
                        #temp_df_1 = pd.DataFrame(list(df.loc[[term]].values))
                        #same_df = same_df.append(temp_df_1)
            diff_file.close()
            same_file.close()
            print("{} finished".format(name))  
        else:
            print("Processing {}".format(name))
            diff_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name" + "/" + name, "w+")
            same_file = open(Path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_same_name" + "/" + name, "w+")
            column_names = list(df.columns.values)
            diff_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            same_file.write("{},{},{},{}\n".format(column_names[0],column_names[1],column_names[2],column_names[3]))
            for term in range(len(gene_name_new)):
                if gene_name_new[term] == "NONE":
                    pass 
                elif other_name_new[term] == "NONE":
                    pass   
                else: 
                    decision = same_filter(other_name_new[term],gene_name_new[term])
                    if decision == "No":
                        ID_list = df.loc[[term]].values
                        ID_list = list(ID_list[0])
                        diff_file.write("{},{},{},{}\n".format(ID_list[0],ID_list[1],ID_list[2],ID_list[3]))
                        #temp_df = pd.DataFrame(list(df.loc[[term]].values), list(df.columns.values)) #change column names of dataframe
                        #sub_df = sub_df.append(temp_df, ignore_index=True) #output the different gene names in to csv files
                    else:
                        ID_list_1 = df.loc[[term]].values
                        ID_list_1 = list(ID_list_1[0])
                        same_file.write("{},{},{},{}\n".format(ID_list_1[0],ID_list_1[1],ID_list_1[2],ID_list_1[3]))
                        #temp_df_1 = pd.DataFrame(list(df.loc[[term]].values))
                        #same_df = same_df.append(temp_df_1)
            diff_file.close()
            same_file.close()
            print("{} finished".format(name))  
        #sub_df = pd.DataFrame(sub_df,columns=list(df.columns.values))
        #sub_df = sub_df.rename(index=str, columns={0: df.columns.values[0], 1: df.columns.values[1], 2: df.columns.values[2], 3:df.columns.values[3]}) # rename columns 
        #sub_df = sub_df.drop_duplicates()
        #same_df = same_df.rename(index=str, columns={0: df.columns.values[0], 1: df.columns.values[1], 2: df.columns.values[2], 3:df.columns.values[3]}) # rename columns 
        #same_df = same_df.drop_duplicates() 
        #sub_df.to_csv(path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_different_name" + "/" + name, index = None, header=True)  #output as csv file
        #same_df.to_csv(path2 + "/" + "Proccessed_Data" + "/" + "Orthologues_with_same_name" + "/" + name, index = None, header=True)
        #printProgressBar(i, len(folder), prefix = 'Progress:', suffix = 'Complete', length = 50)


start = time.time()
makedir(Path2)                    

file_orth_l = os.listdir(Path1)
#for f in file_orth_l:
file_length = os.listdir(Path1 + "/" + file_orth_l[0])
pool = mp.Pool(processes = 7)
pool.map(same_gene,range(len(file_length)))
end = time.time()-start
print("Cost {} min".format(end/60))   
