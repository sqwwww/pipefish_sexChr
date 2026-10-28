

removeIDFile_R1=remove_readsName.m.R1.txt
removeIDFile_R2=remove_readsName.m.R2.txt

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/cDNA/Mbtcombined_S1_L001_1.fq.gz cDNA/MbtcombinedFilt_S1_L001_1.fq.gz $removeIDFile_R1

sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/cDNA/Mbtcombined_S1_L001_2.fq.gz cDNA/MbtcombinedFilt_S1_L001_2.fq.gz $removeIDFile_R2


# test
#removeIDFile_R1=test/test.removeID.txt
#sbatch filter_single_fastq.sh /data/home/zhou3lab/suqianwei/project/seahorse/01.data/10.singleCell/Mbr_testis/0402_1other/combined_data/cDNA/Mbtcombined_S1_L001_1.fq.gz cDNA/MbtcombinedFilt_S1_L001_1.fq.gz $removeIDFile_R1
