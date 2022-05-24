/**
 * Author: Kun Sun (sunkun@szbl.ac.cn)
 * Date: Jun, 2020 
 * This program is part of the Ktrim package
**/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>
#include "common.h"

using namespace std;

//static char tmp_r1[ADAPTER_BUFFER_SIZE];
//static char tmp_r2[ADAPTER_BUFFER_SIZE];
static char tmp_index1[ADAPTER_INDEX_SIZE+1];
static char tmp_index2[ADAPTER_INDEX_SIZE+1];
static char tmp_index3[ADAPTER_INDEX_SIZE+1];

// default parameters
void init_param( ktrim_param &kp ) {
	kp.FASTQ1 = NULL;
	kp.FASTQ2 = NULL;
	kp.outpre = NULL;

	kp.thread = 6;
	kp.min_length = 36;
	kp.phred   = 33;
	kp.minqual = 20;
	kp.quality = 53;
	kp.window  = 5;

	kp.seqKit = NULL;	
	kp.seqA = NULL;
	kp.seqB = NULL;

	kp.adapter_r1 = NULL;
	kp.adapter_r2 = NULL;
	kp.adapter_len= 0;
	kp.adapter_index1 = NULL;
	kp.adapter_index2 = NULL;
	kp.adapter_index3 = NULL;
}

// process user-supplied parameters
int process_cmd_param( int argc, char * argv[], ktrim_param &kp ) {
	const char * prg = argv[0];
	int index;
	int ch;
	while( (ch = getopt(argc, argv, param_list) ) != -1 ) {
		switch( ch ) {
			case '1': kp.FASTQ1 = optarg; break;
			case '2': kp.FASTQ2 = optarg; break;
			case 'o': kp.outpre = optarg; break;
			case 't': kp.thread = atoi(optarg); break;
			case 'k': kp.seqKit = optarg; break;
			case 's': kp.min_length = atoi(optarg); break;
			case 'p': kp.phred = atoi(optarg); break;
			case 'q': kp.minqual = atoi(optarg); break;
			case 'w': kp.window  = atoi(optarg); break;
			case 'a': kp.seqA = optarg; break;
			case 'b': kp.seqB = optarg; break;

			case 'h': usage(); return 100;
			case 'v': cout << VERSION << '\n'; return 100;

			default:
				//cerr << "\033[1;31mError: invalid argument ('" << (char)optopt << "')!\033[0m\n";
				usage();
				return 2;
		}
	}

	argc -= optind;
	argv += optind;
	if( argc > 0 ) {
		cerr << "\033[1;31mError: unrecognized parameter ('" << argv[0] << "')!\033[0m\n";
		usage();
		return 2;
	}

	// check compulsory paramaters
	if( kp.FASTQ1==NULL || kp.FASTQ2==NULL) {
		cerr << "\033[1;31mError: either -1 of -2 are not set!\033[0m\n";
		usage();
		return 2;
	}
	if( kp.outpre == NULL ) {
		cerr << "\033[1;31mError: No output file specified ('-o' not set)!\033[0m\n";
		usage();
		return 3;
	}
	
	// check optional parameters
	if( kp.thread == 0 ) {
		kp.thread = 6;
	} else if ( kp.thread < 4 ) {
		cerr << "\033[1;31mError: thread could not be less than 4!\033[0m\n";
		usage();
		return 4;
	}

//	if( kp.thread > 4 ) {
//		cerr << "Warning: we strongly discourage the usage of more than 4-threads.\n";
//	}
	if( kp.min_length < 10 ) {
		cerr << "\033[1;31mError: invalid min_length! Must be a positive number larger than 10!\033[0m\n";
		usage();
		return 11;
	}
	if( kp.phred==0 || kp.minqual==0 || kp.window==0 ) {
		cerr << "\033[1;31mError: invalid phred/quality/window parameter! Must be positive numbers!\033[0m\n";
		usage();
		return 12;
	}
	kp.quality = kp.phred + kp.minqual;

	// settings of adapters
	if( kp.seqKit != NULL ) {	// use built-in adaptors
		if( strcmp(kp.seqKit, "illumina")==0 || strcmp(kp.seqKit, "ILLUMINA")==0 || strcmp(kp.seqKit, "Illumina")==0 ) {
			kp.adapter_r1 = illumina_adapter_r1;
			kp.adapter_r2 = illumina_adapter_r2;
			kp.adapter_len = illumina_adapter_len;
			kp.adapter_index1 = illumina_index1;
			kp.adapter_index2 = illumina_index2;
			kp.adapter_index3 = illumina_index3;
		} else if( strcmp(kp.seqKit, "Nextera") == 0 || strcmp(kp.seqKit, "nextera") == 0 ) {
			kp.adapter_r1 = nextera_adapter_r1;
			kp.adapter_r2 = nextera_adapter_r2;
			kp.adapter_len = nextera_adapter_len;
			kp.adapter_index1 = nextera_index1;
			kp.adapter_index2 = nextera_index2;
			kp.adapter_index3 = nextera_index3;
		} else if( strcmp(kp.seqKit, "Transposase") == 0 || strcmp(kp.seqKit, "transposase") == 0 ) {
			kp.adapter_r1 = transposase_adapter_r1;
			kp.adapter_r2 = transposase_adapter_r2;
			kp.adapter_len = transposase_adapter_len;
			kp.adapter_index1 = transposase_index1;
			kp.adapter_index2 = transposase_index2;
			kp.adapter_index3 = transposase_index3;
		} else if( strcmp(kp.seqKit, "BGI") == 0 || strcmp(kp.seqKit, "bgi") == 0 ) {
			kp.adapter_r1 = bgi_adapter_r1;
			kp.adapter_r2 = bgi_adapter_r2;
			kp.adapter_len = bgi_adapter_len;
			kp.adapter_index1 = bgi_index1;
			kp.adapter_index2 = bgi_index2;
			kp.adapter_index3 = bgi_index3;
		} else {
			cerr << "\033[1;31mError: unacceptable sequencing kit types!\033[0m\n";
			usage();
			return 13;
		}
		if( kp.seqA!=NULL || kp.seqB!=NULL )
			cerr << "\033[1;32mWarning: '-m' is set! '-a'/'-b' will be ignored!\033[0m\n";
	} else {
		if( kp.seqA==NULL ) {	// no adapters set, use default; note that if -b is set BUT -a is not, report an error
			if( kp.seqB != NULL ) {
				cerr << "\033[1;31mError: '-b' is set while '-a' is not set!\033[0m\n";
				usage();
				return 14;
			}
			kp.adapter_r1 = illumina_adapter_r1;
			kp.adapter_r2 = illumina_adapter_r2;
			kp.adapter_len = illumina_adapter_len;
			kp.adapter_index1 = illumina_index1;
			kp.adapter_index2 = illumina_index2;
			kp.adapter_index3 = illumina_index3;
		} else {
			if( kp.seqB == NULL ) { // seqB not set, then set it to seqA
				kp.seqB = kp.seqA;
				if( kp.FASTQ2 != NULL ) // paired end
					cerr << "\033[1;32mWarning: '-b' is not set for Paired-End data! I will use '-a' as adapter for read 2!\033[0m\n";
			}
			register unsigned int i = 0;
			while( true ) {
				if( kp.seqA[i]=='\0' || kp.seqB[i]=='\0' )
					break;

				++ i;
			}
			if( i < MIN_ADAPTER_SIZE || i > MAX_ADAPTER_SIZE ) {
				cerr << "\033[1;31mError: adapter size must be between " << MIN_ADAPTER_SIZE << " to "
					 << MAX_ADAPTER_SIZE << " bp!\033[0m\n";
				return 20;
			}
			kp.adapter_len = i;

			for( i=0; i!=ADAPTER_INDEX_SIZE; ++i ) {
				tmp_index1[i] = kp.seqA[i];
				tmp_index2[i] = kp.seqB[i];
				tmp_index3[i] = kp.seqA[i+ADAPTER_INDEX_SIZE];
			}
			tmp_index1[ADAPTER_INDEX_SIZE] = '\0';
			tmp_index2[ADAPTER_INDEX_SIZE] = '\0';
			tmp_index3[ADAPTER_INDEX_SIZE] = '\0';

			kp.adapter_r1 = kp.seqA;
			kp.adapter_r2 = kp.seqB;
			kp.adapter_index1 = tmp_index1;
			kp.adapter_index2 = tmp_index2;
			kp.adapter_index3 = tmp_index3;
		}
	}

	return 0;
}

void usage() {
cerr << "\n\033[1;34mUsage: Ktrim [options] -1 Read1.fq -2 Read2.fq -o out.prefix\033[0m\n\n"

	 << "This is a MODIFIED version of Ktrim for Microcket, and is NOT supposed to be called manually by the users.\n\n"

	 << "Compulsory parameters:\n\n"

	 << "  -1 R1.fq[.gz]     Specify the path to the files containing read 1\n"
	 << "  -2 R2.fq[.gz]     Specify the path to the files containing read 2\n\n"
	 << "                    If you have multiple files for your sample, use '"
								<< FILE_SEPARATOR << "' to separate them,\n"
	 << "                    and make sure that all the files are well paired in '-1' and '-2' options\n\n"
     << "                    Note that Gzip-compressed files are supported from version 1.2\n\n"

	 << "  -o out.prefix     Specify the prefix of the output files\n"
	 << "                    Note that output files include trimmed reads in FASTQ format and statistics\n\n"

	 << "Optional parameters:\n\n"


	 << "  -t threads        Specify how many threads should be used (default: 6, must be >=4)\n"
	 << "                    You can set '-t' to 0 to use all threads (automatically detected)\n\n"

	 << "  -p phred-base     Specify the baseline of the phred score (default: 33)\n"
	 << "  -q score          The minimum quality score to keep the cycle (default: 20)\n"
	 << "                    Note that 20 means 1% error rate, 30 means 0.1% error rate in Phred\n\n"
	 << "  -w window         Set the window size for quality check (default: 5)\n"
	 << "                    Ktrim stops when all bases in a consecutive window pass the quality threshold\n\n"

	 << "  -s size           Minimum read size to be kept after trimming (default: 36; must be larger than 10)\n\n"

	 << "  -k kit            Specify the sequencing kit to use built-in adapters\n"
	 << "                    Currently supports 'Illumina' (default), 'Nextera', and 'BGI'\n"
	 << "  -a sequence       Specify the adapter sequence in read 1\n"
	 << "  -b sequence       Specify the adapter sequence in read 2\n"
	 << "                    If '-a' is set while '-b' is not, I will assume that read 1 and 2 use same adapter\n"
	 << "                    Note that '-k' option has a higher priority (when set, '-a'/'-b' will be ignored)\n\n"

	 << "  -h                Show the help information and quit (exit code=0)\n"
	 << "  -v                Show the software version and quit (exit code=0)\n\n"

	 << "Please refer to README.md file for more information (e.g., setting adapters).\n\n"
	 << "\033[1;34mKtrim: extra-fast and accurate adapter- and quality-trimmer.\033[0m\n\n";
}

