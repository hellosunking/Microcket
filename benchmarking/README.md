# Scripts used during benchmark evaluations

Most scripts require parameters; the users could call them without parameter to see the usage.

## The following scripts are used to run the software
- run.juicer.sh
- run.HiC-Pro.sh
- run.HiCUP.sh
- run.fanc.sh
- run.Distiller.sh

Files that were used to count running time are indicated in these scripts.

### IMPORTANT NOTE: the users need to manually amend the paths to the software, genome indices, etc., to run these scripts.

## The following scripts are used to extract the statistics
- make.microcket.stat.pl
- make.juicer.stat.pl
- make.HiC-Pro.stat.pl
- make.HiCUP.stat.pl
- make.fanc.stat.pl
- make.Distiller.stat.pl

## The following script is used to check the consistency and distributions in 2 "pairs" files
- check.consistency.pl (Please note that this program requires LOTS of memory)

## The following scripts are used to call loops and draw venn-plots on IMR90 cell line Hi-C data
- loop/call.loop.sh
- loop/plot.loop.overlap.R

## The following scripts are used to perform *in silico* simulation using [sim3C](https://github.com/cerebis/sim3C "sim3C")
- simulation/generate.reads.sh
- simulation/split.sim3C.pl
- simulation/check.accuracy.pl

## Other scripts and files
- filter.Juicer.pl: used to filter Juicer's result and generate a ".pairs" file
- hicup2pairs: convert HiCUP's output (`sam` format) to `pairs` format (uses LOTS of memory)
- example.fq.list: used by `Microcket`, `run.juicer.sh`, `run.fanc.sh`, and `run.Distiller.sh`
- example.HiC-Pro.conf: used by `run.HiC-Pro.sh`
- example.hicup.conf: used by `run.HiCUP.sh`
- extract.running.time.pl: calculate the time interval among various files
- check.mem.sh: check the memory usage of a program periodically
