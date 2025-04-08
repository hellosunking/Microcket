# Microcket: an extra-fast and vesatile tool for analysis of 3D genomics data (Hi-C, Micro-C, ChIA-PET, and derivant protocols)
Version 1.4, Apr 2025<br />
Authors: Yu Zhao, Mengqi Yang, Qin Peng, Leina Lu, Xiaowen Lyu, and Kun Sun \(sunkun@szbl.ac.cn\)<br />
<br />
Distributed under the
[GNU General Public License v3.0 \(GPLv3\)](https://www.gnu.org/licenses/gpl-3.0.en.html "GPLv3")
for personal and academic usage only.<br />
For detailed information please refer to the license files under `license` directory.

---

## Installation
`Microcket` is written mainly in `C++` for GNU Linux/Unix platforms. After uncompressing the source package, installation
is finished. `Microcket` does not require any specific hardware or OS. The current version has been tested on CentOS v7.5
with Linux kernel v3.10.0.

We recommand the users download the release packages:
```
wget https://github.com/hellosunking/Microcket/archive/refs/tags/v1.4.tar.gz
tar zxf Microcket-1.4.tar.gz
```
Now you will get a new directory named `Microcket-1.3.1`. Most of the required files are included, but you need to build
the genome indices (see the following sections) before you can use `Microckit`.

`Microcket` depends on the following tools:
- [Ktrim](https://github.com/hellosunking/Ktrim "Ktrim")
- [FLASH](http://ccb.jhu.edu/software/FLASH/ "FLASH")
- [BWA](https://github.com/lh3/bwa "BWA")
- [Samtools](http://www.htslib.org/ "Samtools")
- [JuicerTools](https://github.com/aidenlab/JuicerTools "JuicerTools")

The following tools are optional:
- [STAR](https://github.com/alexdobin/STAR "STAR")
- [Pairix](https://github.com/4dn-dcic/pairix "Pairix")
- [cooler](https://github.com/open2c/cooler "cooler")

Pre-compiled executables for these tools are packaged in `bin/` directory (compiled with `g++ v4.8.5` and linked
with `libz v1.2.7`. If you could not run them (which is usually caused by low version of `libc++` or `libz` library),
you may re-compile these programs yourself and replace the ones in the `bin` directory (recommended), then re-compile
the in-house progams for `Microcket` via:
```
user@linux$ make clean && make
```

Note that if you want to generate `.cool` format results, you need to install the
[cooler](https://github.com/open2c/cooler "cooler") package and make sure that it can be called directly from the
command line (i.e., `which cooler` command returns its path).

## Genome indices
Before run `Microcket`, genome index must be built and several annotation files are also needed.
(Change to `Microcket-1.4` directory before running the following commands).

To build genome indices, you need to prepare the genome sequences in fasta format. For instance, if you want to build
indices for human genome hg38, you can download it from the UCSC genome browser:
```
wget -O hg38.p13.fa.gz https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/p13/hg38.p13.fa.gz
gzip -d hg38.p13.fa.gz
```
Since most reference genomes contain unplaced contigs, e.g., hs8_KZ208915v1_fix, we provide a script to remove such
contigs if you do not want them during the analysis (optional):
```
perl util/clean.genome.pl hg38.p13.fa > hg38.p13.clean.fa
```
To build the index for `BWA` (if you want to use `BWA` for the analysis):
```
bin/bwa index -a bwtsw -p index/hg38/BWA/hg38 hg38.p13.fa
## or
## bin/bwa index -a bwtsw -p index/hg38/BWA/hg38 hg38.p13.clean.fa
```
To build the index for `STAR` (if you want to use `STAR` for the analysis):
```
bin/STAR --runThreadN 16 --runMode genomeGenerate --genomeDir index/hg38/STAR --genomeFastaFiles hg38.p13.fa
```

Besides genome indices, you also need to prepare `XXX.info` and `XXX.sam.header` files and put them under the `anno`
directory (`XXX` is the genome id). `XXX.info` is a 2-column file recording the lengths of each chromosome in the genome,
and `XXX.sam.header` is header for SAM/BAM format. Please note that `Microcket` has already packaged such files for `hg38`
and `mm10` therefore you do not need to do this if you are working on these species, and you can refer them as templates
for other species/genomes.<br />

Moreover, we have prepared a utility program `util/build.index.sh` for this task. You need to prepare a genome sequence
in fasta format and run this program to build indices:
```
sh util/build.index.sh <GENOME.FA> <GENOME.ID>
```
The first parameter is the path to the genome sequence file and the second paramter is the identifier of this genome that
you want to use. Note that this script would build indices for both BWA and STAR by default, and you can specify the aligner
as the 3rd parameter after GENOME.ID is you only want to build the indices for only one of them.

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

Authors : Yu Zhao, Mengqi Yang, Fanglei Gong, Qin Peng, Leina Lu, Xiaowen Lyu, and Kun Sun
Software: Kun Sun (sunkun@szbl.ac.cn)
Version : 1.4, Apr 2025

Microcket is an extra-fast and flexible toolkit for Hi-C/Micro-C data analysis.
It has been specifically optimized for long-cycle (100 or longer) Micro-C data.

Parameters:
  -i fq.list   Specify an input file recording paths to fastq files
               Multiple lanes of data are supported; gzip-ed files are also supported.

  -o sid       Set sample id, which will be used as prefix of all output files.

Options:
  -b           Set this flag to keep inter-lane duplications (e.g., biological replicates)
  -k kit       Set sequencing kit for adaptor-trimming. Default: illumina.
                   Supports bgi, illumina, nextera
  -m mode      Set read stitch mode. Default: auto.
                   Supports yes, no, auto

  -g genome    Set genome. Default hg38.
  -a aligner   Set underline aligner. Default: BWA
                   Supports STAR, BWA, STAR-BWA, BWA-STAR
  -e VALUE     Set the minimum alignable proportion of a read to consider as complete mapping
                   Default: 0.5. Must be a number between 0 and 1.
  -Q VALUE     Set the minimum mapping quality for a read to be considered. Default: 10

  -r r1[,r2]   Set resolutions in hic file (use ',' to separate multiple values)
                   Default: 2500000,1000000,500000,250000,100000,50000,25000,10000,5000
                   If the parameter starts with ',', it will be added to the default resolutions
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
reads while BWA for the unstitched ones, and "BWA-STAR" means the opposite manner. In most scenarios, BWA
is recommended as the aligner; if the users want a quick examination of their data (e.g., from shallow
sequencing reads to evaluate the quality of the library), STAR might be used.

The input of `Microcket` is a file containing the paths to the paired-end fastq files, where you can add
as many lines as you like (e.g., 1 line for 1 lane). If your data is generated from biological replicates,
you can set '-b' option to preserve identical reads between each replicates.

### Example 1
Your data is generated using Illumina platform and you have 2 paired-end files (e.g., sequenced twice),
then you can prepare a `fq.list.example1` file as follows:

```
# lines starting with # are considered as comments and ignored
# for each line, only the first 2 columns will be used
# column 1 should record the path to read1
# column 2 should record the path to read2 (paired to read1)
# absolute paths are recommended for the fastq files
# you ca use ',' to list more than 1 files if needed
/path/to/seq1.read1.fq.gz	/path/to/seq1.read2.fq.gz
/path/to/seq2.read1.fq.gz	/path/to/seq2.read2.fq.gz
```

Suppose your data is for human, and your sample id is `test.sample1`, then you can run `Microcket` using
the following command:

```
user@linux$ microcket -i /path/to/fq.list.example1 -o test.sample1
```

### Example 2
You are working on mouse cells and you want to use `mm10` as the reference genome; you have constructed 3
biological replicates and sequenced them on a BGI sequencer (and replicate 1 has 2 paired-end files);
you want to use `STAR` as the underline aligner; you want to visualize the `hic` result in UCSC genome browser
and a `cool` format result for other tools; you want to use 16 CPUs in your computing server; you have prepared
a `fq.list.example2` as follows:

```
# fq.list example 2
/path/to/rep1.lane1.read1.fq.gz,/path/to/rep1.lane2.read1.fq.gz	/path/to/rep1.lane2.read2.fq.gz,/path/to/rep1.lane2.read2.fq.gz
/path/to/rep2.read1.fq.gz	/path/to/rep2.read2.fq.gz
/path/to/rep3.read1.fq.gz	/path/to/rep3.read2.fq.gz
```

then you can run the analysis using the following command:
```
user@linux$ microcket -g mm10 -a star -k bgi -t 16 -buc -i /path/to/fq.list.example2 -o test.sample2
```

## Testing dataset and outputs explanation
As most real HiC/Micro-C datasets are very large, we therefore could not include such data in this source package.
For testing purpose, we prepared a script `util/run.testing.dataset.sh` that will download a small Hi-C dataset containing
~14 million reads from [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/view/SRR4094729 "ENA"), build
hg38 index (if it doesnot exist), and run Microcket automatically. On this dataset, `Microcket` should output the
interaction pairs within 5 minutes and finish the whole analysis (with `hic` and `bam` files generated) in ~10 minutes
using 8 threads on a common computing machine, but the file-downloading time may vary depending on your network speed,
and the index-building step may take ~1 hour using 16 threads.

The following files would be expected after running the script:
```
Microcket.SRR4094729.bwa.unc.log
Microcket.SRR4094729.final.pairs
Microcket.SRR4094729.final.stat
Microcket.SRR4094729.hic
Microcket.SRR4094729.juicer.log
Microcket.SRR4094729.rmdup.log
Microcket.SRR4094729.trim.log
Microcket.SRR4094729.unc2pairs.log
Microcket.SRR4094729.valid.bam
Microcket.SRR4094729.valid.bam.bai
```
In these files, `XXX.valid.bam` (XXX is the parameter in "-o" option) contains the final mappable reads in `BAM` format
(with an index; them bam file will not be generated if '-x' is set), `XXX.final.pairs` records the called `pairs`,
`XXX.hic` records the matrix in `hic` format (optionally in an additional `cool` format), and `XXX.final.stat` records
the key statistics of the analysis (e.g., for QC). For this testing dataset, the statistics would look like this
(note that the numbers may vary when using different versions of Microcket):
```
#Category       Count   Fraction(%)
## Preprocessing and alignment
Total   14,687,360      100.0
Ktrim   12,191,268      83.0
Unique  10,752,452      88.2
Mappable        8,977,487       83.5
## Interactions
Uncalled        171,841 1.9
  Incomplete-mapping    147,614 1.6
  Too-many-segments     0       0.0
  Unpairable    23,265  0.3
  Self-circle   962     0.0
Reported        8,805,646       98.1
  Cis(<1K)      4,371,547       49.6
  Cis(1-10K)    173,137 2.0
  Cis(>=10K)    2,019,269       22.9
  Trans 2,241,693       25.5
```

For comprehensive performance evaluations, we suggest the users use public datasets from literature or consortiums,
e.g., [Rao et al. Cell 2014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525 "Rao et al. Cell 2014"),
[4D nucleome project](https://data.4dnucleome.org "4DN project"), and
[ENCODE project](https://www.encodeproject.org/search/?type=Experiment&assay_title=Hi-C "ENCODE").

---
Please send bug reports to Kun Sun \(sunkun@szbl.ac.cn\).<br />
Microcket is freely available at
[https://github.com/hellosunking/Microcket/](https://github.com/hellosunking/Microcket/ "Microcket @ Github").

