import csv
import os

with open("../am_batch_align_10_25_2021.txt",'r') as name_file:
    f_names = csv.reader(name_file,delimiter='\t')
    for line in f_names:
        #op_file = line[0][0:5]
        #str1 = "(hisat2 -p 4 --dta -x "+line[2]+" -1 "+line[0]+" -2 "+line[1]+" -S sam_output/"+op_file+".sam) >& sam_output/align_"+op_file+"_log.txt &"
        str1 = "(samtools sort -@ 4 -o sorted_bams_10_27_2021/"+line[0]+"_sorted.bam "+line[0]+".sam) >& log_sorting_10_27_2021/log_sort_"+line[0]+" &"
        #str2 = "ln -s "+line[4]+" "+line[5]
        
        print(str1+"\n")
        #print(str2+"\n")

        #os.system (str1)
        #os.system (str2)        
