#!/bin/bash
#
# Author: Kun Sun @ SZBL (sunkun@szbl.ac.cn)
# Date  : Oct 16, 2022
#
set -o nounset
set -o errexit
#command || { echo "command failed"; exit 1; }

echo
echo "This program is designed to run Microcket on a small dataset, for testing purpose only."
echo

PRG=`dirname $0`
mkdir -p testing

## check sequencing data
## SRR4094729: from GSM2297201, Hi-C on STL011.LI11.ADS194
if [ ! -s testing/SRR4094729_1.fastq.gz ]
then
	echo "INFO: sequencing data (SRR4094729_1.fastq.gz) not found, I will download it."
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR409/009/SRR4094729/SRR4094729_1.fastq.gz -o testing/SRR4094729_1.fastq.gz
	if [ $? != 0 ]
	then
		echo ERROR: Download sequencing data failed!
		exit 1
	fi >/dev/stderr
fi

if [ ! -s testing/SRR4094729_2.fastq.gz ]
then
	echo "INFO: sequencing data (SRR4094729_2.fastq.gz) not found, I will download it."
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR409/009/SRR4094729/SRR4094729_2.fastq.gz -o testing/SRR4094729_2.fastq.gz
	if [ $? != 0 ]
	then
		echo ERROR: Download sequencing data failed!
		exit 1
	fi >/dev/stderr
fi
echo -e "$PWD/testing/SRR4094729_1.fastq.gz\t$PWD/testing/SRR4094729_2.fastq.gz" >testing/SRR4094729.fq.list

## check index
if [ ! -s "$PRG/index/STAR/hg38/SA" ]
then
	echo "INFO: index for hg38 is missing, I will try to build one." >/dev/stderr
	if [ ! -s testing/hg38.fa ]
	then
		if [ ! -s testing/hg38.fa.gz ]
		then
			echo "INFO: genome file for hg38 is missing, I will download it." >/dev/stderr
			curl https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -o testing/hg38.fa.gz
			if [ $? != 0 ]
			then
				echo ERROR: Download genome file failed!
				exit 10
			fi >/dev/stderr
		else
			zcat testing/hg38.fa.gz >testing/hg38.fa
		fi
	fi

	$PRG/util/build.index.sh testing/hg38.fa hg38
	if [ $? != 0 ]
	then
		echo ERROR: Build index failed!
		exit 11
	fi > /dev/stderr
fi

## run Microcket
echo "Run Microcket on the testing data ..."
$PRG/microcket -i testing/SRR4094729.fq.list -o testing/Microcket.SRR4094729 -t 8
if [ $? != 0 ]
then
	echo ERROR: Microcket failed!
	exit 20
fi >/dev/stderr

echo "Done: The output files are stored in 'testing/Microcket.SRR4094729'."

