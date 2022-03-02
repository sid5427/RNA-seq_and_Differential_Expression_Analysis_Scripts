import csv
import os


with open("am_batch_align_10_25_2021.txt",'r') as name_file:
    f_names = csv.reader(name_file,delimiter='\t')
    for line in f_names:
        #op_file = line[0][0:5]
        str1 = "(hisat2 -p 4 --dta -x "+line[2]+" -1 "+line[0]+"_1_val_1.fq.gz -2 "+line[0]+"_2_val_2.fq.gz -S sam_output_10_25_2021/"+line[0]+".sam) >& log_sam_output_10_25_2021/align_"+line[0]+"_log.txt &"
        #str1 = "ln -s "+line[2]+" "+line[3]
        #str2 = "ln -s "+line[4]+" "+line[5]
        
        print(str1+"\n")
        #print(str2+"\n")

        #os.system (str1)
        #os.system (str2)        
