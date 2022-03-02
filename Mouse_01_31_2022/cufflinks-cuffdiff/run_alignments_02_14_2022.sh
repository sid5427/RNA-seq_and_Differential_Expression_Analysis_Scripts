#!/bin/bash
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTKO-1_S93_L004_R1_001_val_1.fq.gz -2 MK-CD-PTKO-1_S93_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTKO-1_S93_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTKO-1_S93_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTKO-2_S94_L004_R1_001_val_1.fq.gz -2 MK-CD-PTKO-2_S94_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTKO-2_S94_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTKO-2_S94_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTKO-3_S95_L004_R1_001_val_1.fq.gz -2 MK-CD-PTKO-3_S95_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTKO-3_S95_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTKO-3_S95_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTKO-4_S96_L004_R1_001_val_1.fq.gz -2 MK-CD-PTKO-4_S96_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTKO-4_S96_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTKO-4_S96_L004_log.txt &
##
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTKO-5_S97_L004_R1_001_val_1.fq.gz -2 MK-CD-PTKO-5_S97_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTKO-5_S97_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTKO-5_S97_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-10_S91_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-10_S91_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-10_S91_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-10_S91_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-11_S92_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-11_S92_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-11_S92_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-11_S92_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-1_S82_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-1_S82_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-1_S82_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-1_S82_L004_log.txt
####
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-2_S83_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-2_S83_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-2_S83_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-2_S83_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-3_S84_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-3_S84_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-3_S84_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-3_S84_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-4_S85_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-4_S85_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-4_S85_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-4_S85_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-5_S86_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-5_S86_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-5_S86_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-5_S86_L004_log.txt &
##
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-6_S87_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-6_S87_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-6_S87_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-6_S87_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-7_S88_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-7_S88_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-7_S88_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-7_S88_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-8_S89_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-8_S89_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-8_S89_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-8_S89_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-CD-PTWT-9_S90_L004_R1_001_val_1.fq.gz -2 MK-CD-PTWT-9_S90_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-CD-PTWT-9_S90_L004.sam) >& log_sam_output_02_14_2022/align_MK-CD-PTWT-9_S90_L004_log.txt
####
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-1_S109_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-1_S109_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-1_S109_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-1_S109_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-2_S110_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-2_S110_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-2_S110_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-2_S110_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-3_S111_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-3_S111_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-3_S111_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-3_S111_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-4_S112_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-4_S112_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-4_S112_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-4_S112_L004_log.txt &
##
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-5_S113_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-5_S113_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-5_S113_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-5_S113_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-6_S114_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-6_S114_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-6_S114_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-6_S114_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTKO-7_S115_L004_R1_001_val_1.fq.gz -2 MK-WD-PTKO-7_S115_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTKO-7_S115_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTKO-7_S115_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-10_S107_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-10_S107_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-10_S107_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-10_S107_L004_log.txt
####
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-11_S108_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-11_S108_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-11_S108_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-11_S108_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-1_S98_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-1_S98_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-1_S98_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-1_S98_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-2_S99_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-2_S99_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-2_S99_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-2_S99_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-3_S100_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-3_S100_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-3_S100_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-3_S100_L004_log.txt &
##
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-4_S101_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-4_S101_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-4_S101_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-4_S101_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-5_S102_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-5_S102_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-5_S102_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-5_S102_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-6_S103_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-6_S103_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-6_S103_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-6_S103_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-7_S104_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-7_S104_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-7_S104_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-7_S104_L004_log.txt
####
(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-8_S105_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-8_S105_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-8_S105_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-8_S105_L004_log.txt &

(hisat2 -p 4 --dta-cufflinks -x /home/ssz74/scratch/2022_files/reference_genomes/mouse_GRCm39_01_29_2022/Mus_musculus.GRCm39.dna_sm.primary_assembly -1 MK-WD-PTWT-9_S106_L004_R1_001_val_1.fq.gz -2 MK-WD-PTWT-9_S106_L004_R2_001_val_2.fq.gz -S sam_output_02_14_2022/MK-WD-PTWT-9_S106_L004.sam) >& log_sam_output_02_14_2022/align_MK-WD-PTWT-9_S106_L004_log.txt &
