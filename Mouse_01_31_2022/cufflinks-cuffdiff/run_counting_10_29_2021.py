import csv
import os

with open("../../am_batch_align_10_25_2021.txt",'r') as name_file:
    f_names = csv.reader(name_file,delimiter='\t')
    for line in f_names:
        op_file = line[0][0:5]
        #str1 = "(hisat2 -p 4 --dta -x "+line[2]+" -1 "+line[0]+" -2 "+line[1]+" -S sam_output/"+op_file+".sam) >& sam_output/align_"+op_file+"_log.txt &"
        #str1 = "(samtools sort -@ 4 -o sorted_bams/"+op_file+"_sorted.bam "+line[0]+") >& log_sort_"+op_file+" &"
        str1 = "(htseq-count "+line[0]+"_sorted.bam /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/b73_new_genes_all_mstrg_ids_2.gff3 -i gene_id -t gene -f bam > counts_10_29_2021/"+line[0]+".count.txt) >& logcount_10_29_2021/log_count_"+line[0]+".txt &"        
		#str2 = "ln -s "+line[4]+" "+line[5]
        
        #print(str1+"\n")
        #print(str2+"\n")

        os.system (str1)
        #os.system (str2)        
