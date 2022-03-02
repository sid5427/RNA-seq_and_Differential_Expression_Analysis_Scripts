#!/bin/bash
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABD1_1_val_1.fq.gz -2 ABD1_2_val_2.fq.gz -S sam_output_10_25_2021/ABD1.sam) >& log_sam_output_10_25_2021/align_ABD1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABD2_1_val_1.fq.gz -2 ABD2_2_val_2.fq.gz -S sam_output_10_25_2021/ABD2.sam) >& log_sam_output_10_25_2021/align_ABD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABD3_1_val_1.fq.gz -2 ABD3_2_val_2.fq.gz -S sam_output_10_25_2021/ABD3.sam) >& log_sam_output_10_25_2021/align_ABD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABD4_1_val_1.fq.gz -2 ABD4_2_val_2.fq.gz -S sam_output_10_25_2021/ABD4.sam) >& log_sam_output_10_25_2021/align_ABD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABD5_1_val_1.fq.gz -2 ABD5_2_val_2.fq.gz -S sam_output_10_25_2021/ABD5.sam) >& log_sam_output_10_25_2021/align_ABD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABM1_1_val_1.fq.gz -2 ABM1_2_val_2.fq.gz -S sam_output_10_25_2021/ABM1.sam) >& log_sam_output_10_25_2021/align_ABM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABM2_1_val_1.fq.gz -2 ABM2_2_val_2.fq.gz -S sam_output_10_25_2021/ABM2.sam) >& log_sam_output_10_25_2021/align_ABM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABM3_1_val_1.fq.gz -2 ABM3_2_val_2.fq.gz -S sam_output_10_25_2021/ABM3.sam) >& log_sam_output_10_25_2021/align_ABM3_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABM4_1_val_1.fq.gz -2 ABM4_2_val_2.fq.gz -S sam_output_10_25_2021/ABM4.sam) >& log_sam_output_10_25_2021/align_ABM4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABM5_1_val_1.fq.gz -2 ABM5_2_val_2.fq.gz -S sam_output_10_25_2021/ABM5.sam) >& log_sam_output_10_25_2021/align_ABM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABW1_1_val_1.fq.gz -2 ABW1_2_val_2.fq.gz -S sam_output_10_25_2021/ABW1.sam) >& log_sam_output_10_25_2021/align_ABW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABW2_1_val_1.fq.gz -2 ABW2_2_val_2.fq.gz -S sam_output_10_25_2021/ABW2.sam) >& log_sam_output_10_25_2021/align_ABW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABW3_1_val_1.fq.gz -2 ABW3_2_val_2.fq.gz -S sam_output_10_25_2021/ABW3.sam) >& log_sam_output_10_25_2021/align_ABW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABW4_1_val_1.fq.gz -2 ABW4_2_val_2.fq.gz -S sam_output_10_25_2021/ABW4.sam) >& log_sam_output_10_25_2021/align_ABW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 ABW5_1_val_1.fq.gz -2 ABW5_2_val_2.fq.gz -S sam_output_10_25_2021/ABW5.sam) >& log_sam_output_10_25_2021/align_ABW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFD1_1_val_1.fq.gz -2 AFD1_2_val_2.fq.gz -S sam_output_10_25_2021/AFD1.sam) >& log_sam_output_10_25_2021/align_AFD1_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFD2_1_val_1.fq.gz -2 AFD2_2_val_2.fq.gz -S sam_output_10_25_2021/AFD2.sam) >& log_sam_output_10_25_2021/align_AFD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFD3_1_val_1.fq.gz -2 AFD3_2_val_2.fq.gz -S sam_output_10_25_2021/AFD3.sam) >& log_sam_output_10_25_2021/align_AFD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFD4_1_val_1.fq.gz -2 AFD4_2_val_2.fq.gz -S sam_output_10_25_2021/AFD4.sam) >& log_sam_output_10_25_2021/align_AFD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFD5_1_val_1.fq.gz -2 AFD5_2_val_2.fq.gz -S sam_output_10_25_2021/AFD5.sam) >& log_sam_output_10_25_2021/align_AFD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFM1_1_val_1.fq.gz -2 AFM1_2_val_2.fq.gz -S sam_output_10_25_2021/AFM1.sam) >& log_sam_output_10_25_2021/align_AFM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFM2_1_val_1.fq.gz -2 AFM2_2_val_2.fq.gz -S sam_output_10_25_2021/AFM2.sam) >& log_sam_output_10_25_2021/align_AFM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFM3_1_val_1.fq.gz -2 AFM3_2_val_2.fq.gz -S sam_output_10_25_2021/AFM3.sam) >& log_sam_output_10_25_2021/align_AFM3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFM4_1_val_1.fq.gz -2 AFM4_2_val_2.fq.gz -S sam_output_10_25_2021/AFM4.sam) >& log_sam_output_10_25_2021/align_AFM4_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFM5_1_val_1.fq.gz -2 AFM5_2_val_2.fq.gz -S sam_output_10_25_2021/AFM5.sam) >& log_sam_output_10_25_2021/align_AFM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFW1_1_val_1.fq.gz -2 AFW1_2_val_2.fq.gz -S sam_output_10_25_2021/AFW1.sam) >& log_sam_output_10_25_2021/align_AFW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFW2_1_val_1.fq.gz -2 AFW2_2_val_2.fq.gz -S sam_output_10_25_2021/AFW2.sam) >& log_sam_output_10_25_2021/align_AFW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFW3_1_val_1.fq.gz -2 AFW3_2_val_2.fq.gz -S sam_output_10_25_2021/AFW3.sam) >& log_sam_output_10_25_2021/align_AFW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFW4_1_val_1.fq.gz -2 AFW4_2_val_2.fq.gz -S sam_output_10_25_2021/AFW4.sam) >& log_sam_output_10_25_2021/align_AFW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 AFW5_1_val_1.fq.gz -2 AFW5_2_val_2.fq.gz -S sam_output_10_25_2021/AFW5.sam) >& log_sam_output_10_25_2021/align_AFW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBD1_1_val_1.fq.gz -2 BBD1_2_val_2.fq.gz -S sam_output_10_25_2021/BBD1.sam) >& log_sam_output_10_25_2021/align_BBD1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBD2_1_val_1.fq.gz -2 BBD2_2_val_2.fq.gz -S sam_output_10_25_2021/BBD2.sam) >& log_sam_output_10_25_2021/align_BBD2_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBD3_1_val_1.fq.gz -2 BBD3_2_val_2.fq.gz -S sam_output_10_25_2021/BBD3.sam) >& log_sam_output_10_25_2021/align_BBD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBD4_1_val_1.fq.gz -2 BBD4_2_val_2.fq.gz -S sam_output_10_25_2021/BBD4.sam) >& log_sam_output_10_25_2021/align_BBD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBD5_1_val_1.fq.gz -2 BBD5_2_val_2.fq.gz -S sam_output_10_25_2021/BBD5.sam) >& log_sam_output_10_25_2021/align_BBD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBM1_1_val_1.fq.gz -2 BBM1_2_val_2.fq.gz -S sam_output_10_25_2021/BBM1.sam) >& log_sam_output_10_25_2021/align_BBM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBM2_1_val_1.fq.gz -2 BBM2_2_val_2.fq.gz -S sam_output_10_25_2021/BBM2.sam) >& log_sam_output_10_25_2021/align_BBM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBM3_1_val_1.fq.gz -2 BBM3_2_val_2.fq.gz -S sam_output_10_25_2021/BBM3.sam) >& log_sam_output_10_25_2021/align_BBM3_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBM4_1_val_1.fq.gz -2 BBM4_2_val_2.fq.gz -S sam_output_10_25_2021/BBM4.sam) >& log_sam_output_10_25_2021/align_BBM4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBM5_1_val_1.fq.gz -2 BBM5_2_val_2.fq.gz -S sam_output_10_25_2021/BBM5.sam) >& log_sam_output_10_25_2021/align_BBM5_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBW1_1_val_1.fq.gz -2 BBW1_2_val_2.fq.gz -S sam_output_10_25_2021/BBW1.sam) >& log_sam_output_10_25_2021/align_BBW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBW2_1_val_1.fq.gz -2 BBW2_2_val_2.fq.gz -S sam_output_10_25_2021/BBW2.sam) >& log_sam_output_10_25_2021/align_BBW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBW3_1_val_1.fq.gz -2 BBW3_2_val_2.fq.gz -S sam_output_10_25_2021/BBW3.sam) >& log_sam_output_10_25_2021/align_BBW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBW4_1_val_1.fq.gz -2 BBW4_2_val_2.fq.gz -S sam_output_10_25_2021/BBW4.sam) >& log_sam_output_10_25_2021/align_BBW4_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 BBW5_1_val_1.fq.gz -2 BBW5_2_val_2.fq.gz -S sam_output_10_25_2021/BBW5.sam) >& log_sam_output_10_25_2021/align_BBW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFD1_1_val_1.fq.gz -2 BFD1_2_val_2.fq.gz -S sam_output_10_25_2021/BFD1.sam) >& log_sam_output_10_25_2021/align_BFD1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFD2_1_val_1.fq.gz -2 BFD2_2_val_2.fq.gz -S sam_output_10_25_2021/BFD2.sam) >& log_sam_output_10_25_2021/align_BFD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFD3_1_val_1.fq.gz -2 BFD3_2_val_2.fq.gz -S sam_output_10_25_2021/BFD3.sam) >& log_sam_output_10_25_2021/align_BFD3_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFD4_1_val_1.fq.gz -2 BFD4_2_val_2.fq.gz -S sam_output_10_25_2021/BFD4.sam) >& log_sam_output_10_25_2021/align_BFD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFD5_1_val_1.fq.gz -2 BFD5_2_val_2.fq.gz -S sam_output_10_25_2021/BFD5.sam) >& log_sam_output_10_25_2021/align_BFD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFM1_1_val_1.fq.gz -2 BFM1_2_val_2.fq.gz -S sam_output_10_25_2021/BFM1.sam) >& log_sam_output_10_25_2021/align_BFM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFM2_1_val_1.fq.gz -2 BFM2_2_val_2.fq.gz -S sam_output_10_25_2021/BFM2.sam) >& log_sam_output_10_25_2021/align_BFM2_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFM3_1_val_1.fq.gz -2 BFM3_2_val_2.fq.gz -S sam_output_10_25_2021/BFM3.sam) >& log_sam_output_10_25_2021/align_BFM3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFM4_1_val_1.fq.gz -2 BFM4_2_val_2.fq.gz -S sam_output_10_25_2021/BFM4.sam) >& log_sam_output_10_25_2021/align_BFM4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFM5_1_val_1.fq.gz -2 BFM5_2_val_2.fq.gz -S sam_output_10_25_2021/BFM5.sam) >& log_sam_output_10_25_2021/align_BFM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFW1_1_val_1.fq.gz -2 BFW1_2_val_2.fq.gz -S sam_output_10_25_2021/BFW1.sam) >& log_sam_output_10_25_2021/align_BFW1_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFW2_1_val_1.fq.gz -2 BFW2_2_val_2.fq.gz -S sam_output_10_25_2021/BFW2.sam) >& log_sam_output_10_25_2021/align_BFW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFW3_1_val_1.fq.gz -2 BFW3_2_val_2.fq.gz -S sam_output_10_25_2021/BFW3.sam) >& log_sam_output_10_25_2021/align_BFW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFW4_1_val_1.fq.gz -2 BFW4_2_val_2.fq.gz -S sam_output_10_25_2021/BFW4.sam) >& log_sam_output_10_25_2021/align_BFW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 BFW5_1_val_1.fq.gz -2 BFW5_2_val_2.fq.gz -S sam_output_10_25_2021/BFW5.sam) >& log_sam_output_10_25_2021/align_BFW5_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBD1_1_val_1.fq.gz -2 CBD1_2_val_2.fq.gz -S sam_output_10_25_2021/CBD1.sam) >& log_sam_output_10_25_2021/align_CBD1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBD2_1_val_1.fq.gz -2 CBD2_2_val_2.fq.gz -S sam_output_10_25_2021/CBD2.sam) >& log_sam_output_10_25_2021/align_CBD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBD3_1_val_1.fq.gz -2 CBD3_2_val_2.fq.gz -S sam_output_10_25_2021/CBD3.sam) >& log_sam_output_10_25_2021/align_CBD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBD4_1_val_1.fq.gz -2 CBD4_2_val_2.fq.gz -S sam_output_10_25_2021/CBD4.sam) >& log_sam_output_10_25_2021/align_CBD4_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBD5_1_val_1.fq.gz -2 CBD5_2_val_2.fq.gz -S sam_output_10_25_2021/CBD5.sam) >& log_sam_output_10_25_2021/align_CBD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBM1_1_val_1.fq.gz -2 CBM1_2_val_2.fq.gz -S sam_output_10_25_2021/CBM1.sam) >& log_sam_output_10_25_2021/align_CBM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBM2_1_val_1.fq.gz -2 CBM2_2_val_2.fq.gz -S sam_output_10_25_2021/CBM2.sam) >& log_sam_output_10_25_2021/align_CBM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBM3_1_val_1.fq.gz -2 CBM3_2_val_2.fq.gz -S sam_output_10_25_2021/CBM3.sam) >& log_sam_output_10_25_2021/align_CBM3_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBM4_1_val_1.fq.gz -2 CBM4_2_val_2.fq.gz -S sam_output_10_25_2021/CBM4.sam) >& log_sam_output_10_25_2021/align_CBM4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBM5_1_val_1.fq.gz -2 CBM5_2_val_2.fq.gz -S sam_output_10_25_2021/CBM5.sam) >& log_sam_output_10_25_2021/align_CBM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBW1_1_val_1.fq.gz -2 CBW1_2_val_2.fq.gz -S sam_output_10_25_2021/CBW1.sam) >& log_sam_output_10_25_2021/align_CBW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBW2_1_val_1.fq.gz -2 CBW2_2_val_2.fq.gz -S sam_output_10_25_2021/CBW2.sam) >& log_sam_output_10_25_2021/align_CBW2_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBW3_1_val_1.fq.gz -2 CBW3_2_val_2.fq.gz -S sam_output_10_25_2021/CBW3.sam) >& log_sam_output_10_25_2021/align_CBW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBW4_1_val_1.fq.gz -2 CBW4_2_val_2.fq.gz -S sam_output_10_25_2021/CBW4.sam) >& log_sam_output_10_25_2021/align_CBW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 CBW5_1_val_1.fq.gz -2 CBW5_2_val_2.fq.gz -S sam_output_10_25_2021/CBW5.sam) >& log_sam_output_10_25_2021/align_CBW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFD1_1_val_1.fq.gz -2 CFD1_2_val_2.fq.gz -S sam_output_10_25_2021/CFD1.sam) >& log_sam_output_10_25_2021/align_CFD1_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFD2_1_val_1.fq.gz -2 CFD2_2_val_2.fq.gz -S sam_output_10_25_2021/CFD2.sam) >& log_sam_output_10_25_2021/align_CFD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFD3_1_val_1.fq.gz -2 CFD3_2_val_2.fq.gz -S sam_output_10_25_2021/CFD3.sam) >& log_sam_output_10_25_2021/align_CFD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFD4_1_val_1.fq.gz -2 CFD4_2_val_2.fq.gz -S sam_output_10_25_2021/CFD4.sam) >& log_sam_output_10_25_2021/align_CFD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFD5_1_val_1.fq.gz -2 CFD5_2_val_2.fq.gz -S sam_output_10_25_2021/CFD5.sam) >& log_sam_output_10_25_2021/align_CFD5_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFM1_1_val_1.fq.gz -2 CFM1_2_val_2.fq.gz -S sam_output_10_25_2021/CFM1.sam) >& log_sam_output_10_25_2021/align_CFM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFM2_1_val_1.fq.gz -2 CFM2_2_val_2.fq.gz -S sam_output_10_25_2021/CFM2.sam) >& log_sam_output_10_25_2021/align_CFM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFM3_1_val_1.fq.gz -2 CFM3_2_val_2.fq.gz -S sam_output_10_25_2021/CFM3.sam) >& log_sam_output_10_25_2021/align_CFM3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFM4_1_val_1.fq.gz -2 CFM4_2_val_2.fq.gz -S sam_output_10_25_2021/CFM4.sam) >& log_sam_output_10_25_2021/align_CFM4_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFW1_1_val_1.fq.gz -2 CFW1_2_val_2.fq.gz -S sam_output_10_25_2021/CFW1.sam) >& log_sam_output_10_25_2021/align_CFW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFW2_1_val_1.fq.gz -2 CFW2_2_val_2.fq.gz -S sam_output_10_25_2021/CFW2.sam) >& log_sam_output_10_25_2021/align_CFW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFW3_1_val_1.fq.gz -2 CFW3_2_val_2.fq.gz -S sam_output_10_25_2021/CFW3.sam) >& log_sam_output_10_25_2021/align_CFW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFW4_1_val_1.fq.gz -2 CFW4_2_val_2.fq.gz -S sam_output_10_25_2021/CFW4.sam) >& log_sam_output_10_25_2021/align_CFW4_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 CFW5_1_val_1.fq.gz -2 CFW5_2_val_2.fq.gz -S sam_output_10_25_2021/CFW5.sam) >& log_sam_output_10_25_2021/align_CFW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBD1_1_val_1.fq.gz -2 DBD1_2_val_2.fq.gz -S sam_output_10_25_2021/DBD1.sam) >& log_sam_output_10_25_2021/align_DBD1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBD2_1_val_1.fq.gz -2 DBD2_2_val_2.fq.gz -S sam_output_10_25_2021/DBD2.sam) >& log_sam_output_10_25_2021/align_DBD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBD3_1_val_1.fq.gz -2 DBD3_2_val_2.fq.gz -S sam_output_10_25_2021/DBD3.sam) >& log_sam_output_10_25_2021/align_DBD3_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBD4_1_val_1.fq.gz -2 DBD4_2_val_2.fq.gz -S sam_output_10_25_2021/DBD4.sam) >& log_sam_output_10_25_2021/align_DBD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBD5_1_val_1.fq.gz -2 DBD5_2_val_2.fq.gz -S sam_output_10_25_2021/DBD5.sam) >& log_sam_output_10_25_2021/align_DBD5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBM1_1_val_1.fq.gz -2 DBM1_2_val_2.fq.gz -S sam_output_10_25_2021/DBM1.sam) >& log_sam_output_10_25_2021/align_DBM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBM3_1_val_1.fq.gz -2 DBM3_2_val_2.fq.gz -S sam_output_10_25_2021/DBM3.sam) >& log_sam_output_10_25_2021/align_DBM3_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBM4_1_val_1.fq.gz -2 DBM4_2_val_2.fq.gz -S sam_output_10_25_2021/DBM4.sam) >& log_sam_output_10_25_2021/align_DBM4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBM5_1_val_1.fq.gz -2 DBM5_2_val_2.fq.gz -S sam_output_10_25_2021/DBM5.sam) >& log_sam_output_10_25_2021/align_DBM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBW1_1_val_1.fq.gz -2 DBW1_2_val_2.fq.gz -S sam_output_10_25_2021/DBW1.sam) >& log_sam_output_10_25_2021/align_DBW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBW2_1_val_1.fq.gz -2 DBW2_2_val_2.fq.gz -S sam_output_10_25_2021/DBW2.sam) >& log_sam_output_10_25_2021/align_DBW2_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBW3_1_val_1.fq.gz -2 DBW3_2_val_2.fq.gz -S sam_output_10_25_2021/DBW3.sam) >& log_sam_output_10_25_2021/align_DBW3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBW4_1_val_1.fq.gz -2 DBW4_2_val_2.fq.gz -S sam_output_10_25_2021/DBW4.sam) >& log_sam_output_10_25_2021/align_DBW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/ref_genome/b73_v4_lastest/Zea_mays.B73_RefGen_v4.dna.toplevel -1 DBW5_1_val_1.fq.gz -2 DBW5_2_val_2.fq.gz -S sam_output_10_25_2021/DBW5.sam) >& log_sam_output_10_25_2021/align_DBW5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFD1_1_val_1.fq.gz -2 DFD1_2_val_2.fq.gz -S sam_output_10_25_2021/DFD1.sam) >& log_sam_output_10_25_2021/align_DFD1_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFD2_1_val_1.fq.gz -2 DFD2_2_val_2.fq.gz -S sam_output_10_25_2021/DFD2.sam) >& log_sam_output_10_25_2021/align_DFD2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFD3_1_val_1.fq.gz -2 DFD3_2_val_2.fq.gz -S sam_output_10_25_2021/DFD3.sam) >& log_sam_output_10_25_2021/align_DFD3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFD4_1_val_1.fq.gz -2 DFD4_2_val_2.fq.gz -S sam_output_10_25_2021/DFD4.sam) >& log_sam_output_10_25_2021/align_DFD4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFD5_1_val_1.fq.gz -2 DFD5_2_val_2.fq.gz -S sam_output_10_25_2021/DFD5.sam) >& log_sam_output_10_25_2021/align_DFD5_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFM1_1_val_1.fq.gz -2 DFM1_2_val_2.fq.gz -S sam_output_10_25_2021/DFM1.sam) >& log_sam_output_10_25_2021/align_DFM1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFM2_1_val_1.fq.gz -2 DFM2_2_val_2.fq.gz -S sam_output_10_25_2021/DFM2.sam) >& log_sam_output_10_25_2021/align_DFM2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFM3_1_val_1.fq.gz -2 DFM3_2_val_2.fq.gz -S sam_output_10_25_2021/DFM3.sam) >& log_sam_output_10_25_2021/align_DFM3_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFM4_1_val_1.fq.gz -2 DFM4_2_val_2.fq.gz -S sam_output_10_25_2021/DFM4.sam) >& log_sam_output_10_25_2021/align_DFM4_log.txt
####
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFM5_1_val_1.fq.gz -2 DFM5_2_val_2.fq.gz -S sam_output_10_25_2021/DFM5.sam) >& log_sam_output_10_25_2021/align_DFM5_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFW1_1_val_1.fq.gz -2 DFW1_2_val_2.fq.gz -S sam_output_10_25_2021/DFW1.sam) >& log_sam_output_10_25_2021/align_DFW1_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFW2_1_val_1.fq.gz -2 DFW2_2_val_2.fq.gz -S sam_output_10_25_2021/DFW2.sam) >& log_sam_output_10_25_2021/align_DFW2_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFW3_1_val_1.fq.gz -2 DFW3_2_val_2.fq.gz -S sam_output_10_25_2021/DFW3.sam) >& log_sam_output_10_25_2021/align_DFW3_log.txt &
#
(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFW4_1_val_1.fq.gz -2 DFW4_2_val_2.fq.gz -S sam_output_10_25_2021/DFW4.sam) >& log_sam_output_10_25_2021/align_DFW4_log.txt &

(hisat2 -p 4 --dta -x /home/ssz74/scratch/full_FR697_trinity/necklace/superTranscriptome/final_rep_set/test_rep_b73_new_align/full_b73_genome+newgenes_mstrg -1 DFW5_1_val_1.fq.gz -2 DFW5_2_val_2.fq.gz -S sam_output_10_25_2021/DFW5.sam) >& log_sam_output_10_25_2021/align_DFW5_log.txt &
####
