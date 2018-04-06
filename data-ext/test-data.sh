#!/bin/bash

cd data
wget https://wupathlabs.wustl.edu/FileShare/VirusSeeker/I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz
tar -xvzf I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz

cd I1164_12629_Harvard_SIV_196_06_2_24_12
zcat I1164_12629_Harvard_SIV_196_06_2_24_12_SE1.fastq.gz | head -8000 | gzip > I1164_12629_Harvard_SIV_196_06_2_24_12_mini_SE1.fastq.gz
zcat I1164_12629_Harvard_SIV_196_06_2_24_12_SE2.fastq.gz | head -8000 | gzip > I1164_12629_Harvard_SIV_196_06_2_24_12_mini_SE2.fastq.gz

mv I1164_12629_Harvard_SIV_196_06_2_24_12_mini_SE{1,2}.fastq.gz ../
cd ..
rm -r I1164_12629_Harvard_SIV_196_06_2_24_12

# rm I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz

