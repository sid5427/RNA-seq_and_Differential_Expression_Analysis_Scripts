# -*- coding: utf-8 -*-
"""
Created on Wed May  9 14:03:49 2018

@author: ssz74
"""

import csv,os

#name_fpkm_dir = "fpkms_files"
#name_fpkm_dir = "cufflinks_07_2019"
name_fpkm_dir = "counts_01_30_2022" #<<---- change to foldername with counts/fpkm files

directory = os.fsencode(name_fpkm_dir)

# List of your files
file_names = []


#files names populated by reading fpkm directory
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #if filename.endswith(".fpkm_tracking.txt"): 
    # print(os.path.join(directory, filename))
    #file_names.append(name_fpkm_dir+"/"+filename)
    file_names.append(name_fpkm_dir+"/"+filename)


#check names by print file names list
print (file_names)
file_names.sort()

# Output list of generator objects
o_data = []

# Open files in the succession and 
# store the file_name as the first
# element followed by the elements of
# the third column.
for afile in file_names:
    file_h = open(afile)
    a_list = []
    a_list.append(afile)
    csv_reader = csv.reader(file_h, delimiter="\t")
    print ("\n"+afile)
    for row in csv_reader:
        a_list.append(row[1])
        
    # Convert the list to a generator object
    o_data.append((n for n in a_list))
    file_h.close()

# Use zip and csv writer to iterate
# through the generator objects and 
# write out to the output file


with open('rna_seq_counts_11_1_2021.txt','w',newline='\n') as op_file:  ##<<--- change to combined matrix file name
    csv_writer = csv.writer(op_file, delimiter='\t')
    for row in list(zip(*o_data)):
        csv_writer.writerow(row)
        
op_file.close()
