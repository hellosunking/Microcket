#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <unordered_set>
#include <ctime>
#include <thread>
using namespace std;
/*
 * Author: Kun Sun @ SZBL
*/

static int write_thread;
const static chrono::milliseconds waiting_time_for_writing(1);
const unsigned int BATCH = 1 << 16;
typedef unsigned long fqKey;

// PE READ
typedef struct {
	string id;
	string seq1;
	string qual1;
	string id2;
	string seq2;
	string qual2;
	bool keep;
} PEREAD;

typedef struct {
	PEREAD *A;
	unsigned int cnt_A;
	PEREAD *C;
	unsigned int cnt_C;
	PEREAD *G;
	unsigned int cnt_G;
	PEREAD *T;
	unsigned int cnt_T;
	unsigned int discard;
} readbatch;

typedef struct {
	char *readx;
	char *outprefix;

	unsigned int hskip1;	// how many to skip in read 1
	unsigned int keylen1;	// the length of the key in read 1
	unsigned int epos1;		// minimum size for read 1

	unsigned int hskip2;	// how many to skip in read 2
	unsigned int keylen2;	// the length of the key in read 2
	unsigned int epos2;		// minimum size for read 2

	FILE *fout1;
	FILE *fout2;
} kparam;

// statistics
typedef struct {
	unsigned int uniq;
	unsigned int dup;
	unsigned int discard;
	char *result1;
	char *result2;
} kstat;

// functions
void usage( const char * prg ) {
	cerr << "\nUsage: " << prg << " [options] -i <interleaved.paired-end.fq> -o <output.prefix>\n"
		 << "\nOptions:\n"
		 << "  -k <int>  Skip the heading cycles in read 1 (default: 5)\n"
		 << "  -K <int>  Skip the heading cycles in read 2 (default: 5)\n"
		 << "  -s <int>  Size of the KEY in read 1 (default: 16)\n"
		 << "  -S <int>  Size of the KEY in read 2 (default: 16)\n"
		 << "\nThis program is designed to remove the duplicate reads from the FASTQ data using 4+1 threads."
		 << "\n\nIMPORTANT NOTEs:"
		 << "\nThe total KEY size in read1 and read2 must >=16 and <=32."
		 << "\nWhen running this program on the adapter-and-quality trimmed data, please mind the read length,"
		 << "\nreads that are shorter than Skip1+Key1 or Skip2+Key2 will be discarded."
		 << "\n\nLog and FASTQ files will be written.\n\n";

	exit(2);
}

int load_batch( ifstream &fq, readbatch *rb, const kparam &kp ) {
	rb->cnt_A=0;
	rb->cnt_C=0;
	rb->cnt_G=0;
	rb->cnt_T=0;
	rb->discard=0;
	int cnt_all = 0;

	string tmp_id, tmp_seq;

	while( true ) {
		getline( fq, tmp_id );
		if( fq.eof() )break;
		getline( fq, tmp_seq );

		char firstLetter = 'N';
		if( tmp_seq.size() >= kp.epos1 ) {
			firstLetter = tmp_seq[ kp.hskip1 ];
		}

		if( firstLetter == 'N' ) {	// discard this one
			++ rb->discard;
			getline( fq, tmp_id );getline( fq, tmp_id );
			getline( fq, tmp_id );getline( fq, tmp_id );getline( fq, tmp_id );getline( fq, tmp_id );
		} else {
			if( firstLetter == 'A' ) {
				rb->A[rb->cnt_A].id = tmp_id;rb->A[rb->cnt_A].seq1 = tmp_seq;
				getline( fq, tmp_id );getline( fq, rb->A[rb->cnt_A].qual1 );
				getline( fq, rb->A[rb->cnt_A].id2 );getline( fq, rb->A[rb->cnt_A].seq2 );
				getline( fq, tmp_id );getline( fq, rb->A[rb->cnt_A].qual2 );
				rb->A[rb->cnt_A].keep = false;
				++ rb->cnt_A;
			} else if ( firstLetter == 'C' ) {
				rb->C[rb->cnt_C].id = tmp_id;rb->C[rb->cnt_C].seq1 = tmp_seq;
				getline( fq, tmp_id );getline( fq, rb->C[rb->cnt_C].qual1 );
				getline( fq, rb->C[rb->cnt_C].id2 );getline( fq, rb->C[rb->cnt_C].seq2 );
				getline( fq, tmp_id );getline( fq, rb->C[rb->cnt_C].qual2 );
				rb->C[rb->cnt_C].keep = false;
				++ rb->cnt_C;
			} else if ( firstLetter == 'G' ) {
				rb->G[rb->cnt_G].id = tmp_id;rb->G[rb->cnt_G].seq1 = tmp_seq;
				getline( fq, tmp_id );getline( fq, rb->G[rb->cnt_G].qual1 );
				getline( fq, rb->G[rb->cnt_G].id2 );getline( fq, rb->G[rb->cnt_G].seq2 );
				getline( fq, tmp_id );getline( fq, rb->G[rb->cnt_G].qual2 );
				rb->G[rb->cnt_G].keep = false;
				++ rb->cnt_G;
			} else {
				rb->T[rb->cnt_T].id = tmp_id;rb->T[rb->cnt_T].seq1 = tmp_seq;
				getline( fq, tmp_id );getline( fq, rb->T[rb->cnt_T].qual1 );
				getline( fq, rb->T[rb->cnt_T].id2 );getline( fq, rb->T[rb->cnt_T].seq2 );
				getline( fq, tmp_id );getline( fq, rb->T[rb->cnt_T].qual2 );
				rb->T[rb->cnt_T].keep = false;
				++ rb->cnt_T;
			}
		}

		++ cnt_all;
		if( cnt_all == BATCH ) break;
	}

	return cnt_all;
}

void do_rmdup( PEREAD *r, int cnt, unordered_set<fqKey> &log, kstat &ks, const kparam &kp, int tn ) {
	const PEREAD *rend = r + cnt;
	int curr_length1 = 0;
	int curr_length2 = 0;

	for( register PEREAD *pread=r; pread!=rend; ++pread ) {
		// check size
		if( pread->seq1.size() < kp.epos1 || pread->seq2.size() < kp.epos2 ) {	// too short
			++ ks.discard;
			continue;
		}

		// make key and check N
		register const char *p = pread->seq1.c_str();
		register const char *q = pread->seq2.c_str();
		bool nflag = false;
		fqKey key = 0;
		for( register int i=kp.hskip1; i!=kp.epos1; ++i ) {
			// make key.r1
			key <<= 2;
			if( p[i] == 'A' || p[i] == 'a' )		key |= 1;
			else if( p[i] == 'T' || p[i] == 't' )	key |= 2;
			else if( p[i] == 'C' || p[i] == 'c' )	key |= 0;
			else if( p[i] == 'G' || p[i] == 'g' )	key |= 3;
			else {	// N encounted
				nflag = true;
				break;
			}
		}
		if( ! nflag ) {
			for( register int i=kp.hskip2; i!=kp.epos2; ++i ) {
				// make key.r2
				key <<= 2;
				if( q[i] == 'A' || q[i] == 'a' )		key |= 1;
				else if( q[i] == 'T' || q[i] == 't' )	key |= 2;
				else if( q[i] == 'C' || q[i] == 'c' )	key |= 0;
				else if( q[i] == 'G' || q[i] == 'g' )	key |= 3;
				else {	// N encounted
					nflag = true;
					break;
				}
			}
		}

		if( nflag ) {
			++ ks.discard;
			continue;
		}

		// now this fragment is a valid one: long enough and do not contain N
		if( log.find( key ) == log.end() ) {	//do not have this key before
			log.emplace( key );
			pread->keep = true;
			++ ks.uniq;

			curr_length1 += sprintf( ks.result1+curr_length1, "%s\n%s\n+\n%s\n",
						pread->id.c_str(),  pread->seq1.c_str(), pread->qual1.c_str() );
			curr_length2 += sprintf( ks.result2+curr_length2, "%s\n%s\n+\n%s\n",
						pread->id2.c_str(), pread->seq2.c_str(), pread->qual2.c_str() );
		} else {
			++ ks.dup;
		}
	}

	// wait for my turn to write output
	while( true ) {
		if( tn == write_thread ) {
			fwrite( ks.result1, 1, curr_length1, kp.fout1 );
			fwrite( ks.result2, 1, curr_length2, kp.fout2 );

			++ write_thread;
			return;
		} else {
			this_thread::sleep_for( waiting_time_for_writing );
		}
	}
}

void init_param( kparam &kp, int argc, char *argv[] ) {
	// init parameters
	kp.hskip1  = 5;	// how many to skip
	kp.keylen1 = 16;// the length of key1
	kp.hskip2  = 5;	// how many to skip
	kp.keylen2 = 16;// the length of key2
//	kp.epos = 21;

	kp.readx = NULL;
	kp.outprefix = NULL;
	string srstr, nrstr, rpstr;
	int opt;
	//char * optvalue;
	while( (opt = getopt(argc, argv, "i:o:k:K:s:S:")) != -1 ) {
		switch( opt ) {
			case 'i': kp.readx=optarg;break;
			case 'o': kp.outprefix=optarg;break;
			case 'k': kp.hskip1=atoi(optarg);break;
			case 'K': kp.hskip2=atoi(optarg);break;
			case 's': kp.keylen1=atoi(optarg);break;
			case 'S': kp.keylen2=atoi(optarg);break;
			case '?':
			default : usage( argv[0] );break;
		}
	}

	// check parameters
	if( kp.readx==NULL || kp.outprefix==NULL ) {
		usage( argv[0] );
	}
	if( kp.keylen1+kp.keylen2>32 || kp.keylen1+kp.keylen2<16) {
		cerr << "Error: invalid key sizes!\n";
		exit(1);
	}
	kp.epos1 = kp.hskip1 + kp.keylen1;
	kp.epos2 = kp.hskip2 + kp.keylen2;

	string outfile = kp.outprefix;
	outfile += ".read1.fq";
	kp.fout1 = fopen( outfile.c_str(), "a" );
	outfile = kp.outprefix;
	outfile += ".read2.fq";
	kp.fout2 = fopen( outfile.c_str(), "a" );
	if( kp.fout1==NULL || kp.fout2==NULL ) {
		cerr << "Error: open output files failed!\n";
		exit(1);
	}
}

int main( int argc, char *argv[] ) {
//	cerr << "Initiate parameters\n";
	kparam kp;
	init_param( kp, argc, argv );

	// prepare files
//	cerr << "Prepare files\n";
	ifstream fqx;
	if( kp.readx[0]=='-' && kp.readx[1]=='\0' ) {
		fqx.open( "/dev/stdin" );
	} else {
		fqx.open( kp.readx );
	}
	if( fqx.fail() ) {
		cerr << "Error: read fastq failed!\n";
		fqx.close();
		return 10;
	}

	// container and statistics
	kstat ks_A, ks_C, ks_G, ks_T;
	ks_A.uniq=0;ks_A.dup=0;ks_A.discard=0;ks_A.result1=new char [ BATCH << 10 ];ks_A.result2=new char [ BATCH << 10 ];
	ks_C.uniq=0;ks_C.dup=0;ks_C.discard=0;ks_C.result1=new char [ BATCH << 10 ];ks_C.result2=new char [ BATCH << 10 ];
	ks_G.uniq=0;ks_G.dup=0;ks_G.discard=0;ks_G.result1=new char [ BATCH << 10 ];ks_G.result2=new char [ BATCH << 10 ];
	ks_T.uniq=0;ks_T.dup=0;ks_T.discard=0;ks_T.result1=new char [ BATCH << 10 ];ks_T.result2=new char [ BATCH << 10 ];

	// prepare memory
//	cerr << "Prepare memory\n";
	readbatch rb1, rb2, *wkr, *ldr;	//A, B is real data, wkr and ldr are indicators for work and load
	rb1.A = new PEREAD [ BATCH ];
	rb1.C = new PEREAD [ BATCH ];
	rb1.G = new PEREAD [ BATCH ];
	rb1.T = new PEREAD [ BATCH ];
	rb2.A = new PEREAD [ BATCH ];
	rb2.C = new PEREAD [ BATCH ];
	rb2.G = new PEREAD [ BATCH ];
	rb2.T = new PEREAD [ BATCH ];

//	cerr << "Loading data\n";
	int nwk = load_batch( fqx, &rb1, kp );
	int discard_in_loading = rb1.discard;
	wkr = &rb1;
	ldr = &rb2;

	unordered_set<fqKey> log_A, log_C, log_G, log_T;

	// start analysis
	bool nextBatch = true;
	unsigned int total = 0;
	int discarded_in_load;
	do {
		write_thread = 0;
		if( fqx.eof() ) {	//no need to load data, use 4 threads to do the analysis
			omp_set_num_threads( 4 );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				if( tn == 0 ) {
					do_rmdup( wkr->A, wkr->cnt_A, log_A, ks_A, kp, tn );
				} else if( tn == 1 ) {
					do_rmdup( wkr->C, wkr->cnt_C, log_C, ks_C, kp, tn );
				} else if( tn == 2 ) {
					do_rmdup( wkr->G, wkr->cnt_G, log_G, ks_G, kp, tn );
				} else {
					do_rmdup( wkr->T, wkr->cnt_T, log_T, ks_T, kp, tn );
				}
			}
			nextBatch = false;
		} else { // split into 2 threads, 1 for loading data, another for analysis
			omp_set_num_threads( 5 );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				if( tn == 0 ) {
					do_rmdup( wkr->A, wkr->cnt_A, log_A, ks_A, kp, tn );
				} else if( tn == 1 ) {
					do_rmdup( wkr->C, wkr->cnt_C, log_C, ks_C, kp, tn );
				} else if( tn == 2 ) {
					do_rmdup( wkr->G, wkr->cnt_G, log_G, ks_G, kp, tn );
				} else if( tn == 3 ) {
					do_rmdup( wkr->T, wkr->cnt_T, log_T, ks_T, kp, tn );
				} else {
					total += load_batch( fqx, ldr, kp );
					discard_in_loading += ldr->discard;
				}
			}
		}
		// update parameters for next loop
		readbatch *tmp = ldr;
		ldr = wkr;
		wkr = tmp;
	} while( nextBatch );
//	cerr << "Done.\n";
	fqx.close();
	fclose( kp.fout1 );
	fclose( kp.fout2 );

	string out = kp.outprefix;
	out += ".log";
	ofstream flog( out.c_str(), ios::app );
	if( flog.fail() ) {
		cerr << "Error: write log failed!\n";
		return 10;
	}

	unsigned int uniq = ks_A.uniq+ks_C.uniq+ks_G.uniq+ks_T.uniq;
	unsigned int dup  = ks_A.dup+ks_C.dup+ks_G.dup+ks_T.dup;
	unsigned int discard = ks_A.discard+ks_C.discard+ks_G.discard+ks_T.discard+discard_in_loading;
	flog << "Total\t"	<< uniq+dup+discard
		 << "\nUniq\t"	<< uniq
		 << "\nDup\t"	<< dup
		 << "\nDiscard\t"<< discard << '\n';
	flog.close();

	delete [] ks_A.result1;delete [] ks_C.result1;delete [] ks_G.result1;delete [] ks_T.result1;
	delete [] ks_A.result2;delete [] ks_C.result2;delete [] ks_G.result2;delete [] ks_T.result2;
	delete [] rb1.A;delete [] rb1.C;delete [] rb1.G;delete [] rb1.T;
	delete [] rb2.A;delete [] rb2.C;delete [] rb2.G;delete [] rb2.T;
	return 0;
}

