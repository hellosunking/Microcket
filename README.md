# Microcket: an extra-fast and flexiable tool for 3D genomics data analysis
Version 1.0.0, May 2022<br />
Authors: Yu Zhao, Mengqi Yang, Qin Peng, Leina Lu, Xiaowen Lyu, and Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please refer to the license files under `license` directory.

---

## Installation
`Microcket` is written mainly in `C++` for GNU Linux/Unix platforms. `Microcket` depends on the following
tools:

- [STAR](https://github.com/alexdobin/STAR "STAR")
- [BWA](https://github.com/lh3/bwa "BWA") (optional)
- [Samtools](http://www.htslib.org/ "Samtools")
- [FLASH](http://ccb.jhu.edu/software/FLASH/ "FLASH")
- [Pairix](https://github.com/4dn-dcic/pairix "Pairix") (optional)
- [JuicerTools](https://github.com/aidenlab/JuicerTools "JuicerTools")

Pre-compiled executables for these tools are included in `bin/` directory (compiled with `g++ v4.8.5` and
linked with `libz v1.2.7`. If you could not run them (which is usually caused by low version of `libc++` or
`libz` library), you may re-compile these programs yourself and replace the ones in the `bin` directory,
then re-compile the in-house progams for `Microcket` via:
```
user@linux$ make clean && make
```

Note that if you want to generate `.cool` format results, you need to install the
[cooler](https://github.com/open2c/cooler "cooler") package and make sure that it can be called directly from
the command line (i.e., `which cooler` command returns its path).

## Pre-requirements
Before run `Microcket`, genome index must be built and several annotation files are also needed.

To build genome indices, you need the genome sequence in fasta format. For instance, if you want to build indices
for human genome hg38, you can download it from the UCSC genome browser:
```
wget -O hg38.p13.fa.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p13/hg38.p13.fa.gz
gzip -d hg38.p13.fa.gz
```
To build the index for `STAR` (if you want to use `STAR` for the analysis):
```
bin/STAR --runThreadN 32 --runMode genomeGenerate --genomeDir index/STAR/hg38 --genomeFastaFiles hg38.p13.fa
```
To build the index for `BWA` (if you want to use `BWA` for the analysis):
```
bin/bwa index -a bwtsw -p index/BWA/hg38 hg38.p13.fa
```
On the other hand, if you already had built such indices before, you can link them to `index/`.<br />

Besides genome indices, you also need to prepare `XXX.info` and `XXX.sam.header` files and put them under the
`anno` directory (`XXX` is the genome id). `XXX.info` is a 2-column file recording the lengths of each chromosome
in the genome, and `XXX.sam.header` is header for SAM/BAM format. Please note that `Microcket` has already
packaged such files for `hg38` and `mm10` therefore you do not need to do this if you are working on these
species, and you can refer them as templates for other species/genomes.

## Run Microcket
The main program is `microcket` under the same directory as this `README.md` file. You can add its path to
your `~/.bashrc` file under the `PATH` variable to call it from anywhere; or you can run the following command
to add it to your current session:
```
user@linux$ PATH=$PATH:$PWD
```

Call `microcket` without any parameters to see the usage (or use '-h' option):
```
Usage: microcket [options] -i <fq.list> -o <sid>

Authors : Yu Zhao, Mengqi Yang, Qin Peng, Leina Lu, Xiaowen Lyu, and Kun Sun
Software: Kun Sun (sunkun@szbl.ac.cn)
Version: 1.0.0, May 2022

Microcket is an extra-fast and flexiable toolkit for Hi-C/Micro-C data analysis.
It has been specifically optimized for long-cycle (100 or longer) Micro-C data.

Options:
  -i fq.list   Specify an input file recording paths to fastq files
               Multiple lanes of data are supported; gzip-ed files are also supported.

  -o sid       Set sample id, which will be used as prefix of all output files.

  -b           Set this flag to keep inter-lane duplications (e.g., biological replicates)
  -m mode      Set read stitch mode. Default: auto.
                   Supports yes, no, auto

  -g genome    Set genome. Default hg38.
  -a aligner   Set underline aligner. Default: STAR.
                   Supports STAR, BWA, STAR-BWA, BWA-STAR
  -k kit       Set sequencing kit for adaptor-trimming. Default: illumina.
                   Supports bgi, illumina, nextera

  -r r1[,r2]   Set resolutions in hic file (use ',' to separate multiple values)
               Default: 2500000,1000000,500000,250000,100000,50000,25000,10000,5000
  -u           Set this flag to generate UCSC-compatible hic files (slower)
  -c           Set this flag to generate cool files (requires cooler)
                   Note: the smallest bin in '-r' will be used during cool file generation.
  -x           Set to skip BAM file generation.

  -t thread    Set running threads.
                   Default: all threads available (>=4 required)

  -q           Set this flag to run in quiet mode (only report fatal errors)
  -v           Show software version and quit
  -h           Show this help information and quit

Microcket is freely available at https://github.com/hellosunking/Microcket
```

Note: Microcket supports 4 modes in alignment in '-a' option; "STAR-BWA" means using STAR for the stitched
reads while BWA for the unstitched ones, and "BWA-STAR" means the opposite manner.

The input of `Microcket` is a file containing the paths to the paired-end fastq files, where you can add
as many lines as you like (e.g., 1 line for 1 lane). If your data is generated from biological replicates,
you can set '-b' option to preserve identical reads between each replicates.

### Example 1
Your data is generated using Illumina platform with 2 lanes, then you can prepare a `fq.list.example1`
file as follows:

```
# lines starting with # are considered as comments
# for each line, only the first 2 columns will be used
# absolute paths are recommended for the fastq files
/path/to/lane1.read1.fq.gz	/path/to/lane1.read2.fq.gz
/path/to/lane2.read1.fq.gz	/path/to/lane2.read2.fq.gz
```

Surppose your data is for human, and your sample id is `test.sample1`, then you can run `Microcket` using
the following command:

```
user@linux$ microcket -i /path/to/fq.list.example1 -o test.sample1
```

### Example 2
You are working on mouse cells and you want to use `mm10` as the reference genome; you have constructed 3
biological replicates and sequenced them on a BGI sequencer (and replicate 1 has two lanes); you want to use
`BWA` as the underline aligner; you want to visualize the `hic` result in UCSC genome browser and a `cool`
format result for other tools; you want to use 16 CPUs in your computing server; you have prepared a
`fq.list.example2` as follows:

```
# fq.list example 2
/path/to/rep1.lane1.read1.fq.gz,/path/to/rep1.lane2.read1.fq.gz	/path/to/rep1.lane2.read2.fq.gz,/path/to/rep1.lane2.read2.fq.gz
/path/to/rep2.read1.fq.gz	/path/to/rep2.read2.fq.gz
/path/to/rep3.read1.fq.gz	/path/to/rep3.read2.fq.gz
```

then you can run the analysis using the following command:
```
user@linux$ microcket -g mm10 -a bwa -k bgi -t 16 -buc -i /path/to/fq.list.example2 -o test.sample2
```

## Outputs explanation
`Microcket` outputs the final mappable reads in BAM format (with an index) unless '-X' is set, called interactions
in `hic` format (optionally in an additional `cool` format), and key statistics of the analysis.

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Microcket is freely available at
[https://github.com/hellosunking/Microcket/](https://github.com/hellosunking/Microcket/ "Microcket @ Github").

