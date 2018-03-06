#!/bin/bash
# download and unzip virusseeker test dataset
cd data
wget https://wupathlabs.wustl.edu/FileShare/VirusSeeker/I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz
tar -xvzf I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz
rm I1164_12629_Harvard_SIV_196_06_2_24_12_rawData.tgz
