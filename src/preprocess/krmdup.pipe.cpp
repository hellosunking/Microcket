#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <getopt.h>
#include <tr1/unordered_set>

using namespace std;
using namespace std::tr1;

/*
 * Author: Kun Sun @ SZBL
*/

const unsigned int BATCH = 1 << 16;
typedef unsigned long fqKey;

// PE READ
typedef struct {
	string id;
	string seq1;
	string seq2;
	string qual1;
	string qual2;
	bool keep;
} PEREAD;

typedef struct {
	char *readx;
	char *outprefix;
	unsigned int hskip;	// how many to skip
	unsigned int keylen;		// the length of the key, IT IS FIXED IN v3
	unsigned int epos;
	bool keepShort;	// if true, keep short reads; false to discard
	bool keepN;	// true to keep the reads that have Ns; false to discard
} kparam;

// statistics
typedef struct {
	unsigned int uniq;
	unsigned int dup;
	unsigned int cN;
	unsigned int cS;	// countN, countShort
} kstat;

// functions
void usage( const char * prg ) {
	cerr << "\nUsage: " << prg << " [options] -i in.fqx -o output.prefix\n";

	cerr << "\nOptions:\n"
		 << "    -k  Skip the heading cycles that usually have lower quality (default: 5)\n"
		 << "    -s  Handle the short reads (default: discard, accept 'discard' or 'keep')\n"
		 << "    -n  Handle the read with Ns (default: discard, accept 'discard' or 'keep')\n";

	cerr << "\nThis program is designed to remove the duplicate reads from the FASTQ data using 2 threads."
		 << "\n\nIMPORTANT NOTE: input is in 5-line fqx format (id,seq1,qual1,seq2,qual2)."
		 << "\nThe KEY size is fixed to 16+16 to balance sensitivity and error-tolerance."
		 << "\nWhen running this program on the adapter-and-quality trimmed data, please mind the read length,"
		 << "\nthe reads that are shorter than skip + len will not be considered, discard or keep them depends on your choice."
		 << "\n\nLog file will be written, while reads wwill be output to STDOUT in interleaved-fastq for FLASH.\n\n";

	exit(2);	// it seems that other people use code 2 as usage error
}

int load_batch( ifstream &fqx, PEREAD *r ) {
	int cnt = 0;
	string tmp;
	while( 1 ) {
		getline( fqx, r[cnt].id    );
		if( fqx.eof() )break;
		getline( fqx, r[cnt].seq1  );
		getline( fqx, r[cnt].qual1 );
		getline( fqx, r[cnt].seq2  );
		getline( fqx, r[cnt].qual2 );
		r[cnt].keep = false;

		++ cnt;
		if( cnt==BATCH ) {
			break;
		}
	}

	return cnt;
}

void do_rmdup( PEREAD *r, int cnt, unordered_set<fqKey> &log, kstat &ks, const kparam &kp ) {
	PEREAD *rend = r + cnt;
	for( PEREAD *pread=r; pread!=rend; ++pread ) {
		// check size
		if( pread->seq1.size() < kp.epos || pread->seq2.size() < kp.epos ) {	// too short
			++ ks.cS;
			if( kp.keepShort ) {
				pread->keep = true;
			}
			continue;
		}

		// make key and check N
		register const char *p = pread->seq1.c_str();
		register const char *q = pread->seq2.c_str();
		bool nflag = false;
		fqKey key = 0;
		for( int i=kp.hskip; i!=kp.epos; ++i ) {
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
			for( int i=kp.hskip; i!=kp.epos; ++i ) {
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
			++ ks.cN;
			if( kp.keepN ) {
				pread->keep = true;
			}
			continue;
		}

		// now this fragment is a valid one: long enough and do not contain N
		if( log.find( key ) == log.end() ) {	//do not have this key before
			log.insert( key );
			pread->keep = true;
			++ ks.uniq;
		} else {
			++ ks.dup;
		}
	}
}

void write_batch( PEREAD *wkr, int nwk ) {
	PEREAD *pend = wkr + nwk;
	for( PEREAD *p=wkr; p!=pend; ++p ) {
		if( p->keep ) {
			cout << p->id << '\n' << p->seq1 << "\n+\n" << p->qual1 << "\n"
				<< p->id << '\n' << p->seq2 << "\n+\n" << p->qual2 << "\n";
		}
	}
}

void init_param( kparam &kp,  int argc, char *argv[] ) {
	// init parameters
	kp.hskip = 5;	// how many to skip
	kp.keylen = 16;		// the length of the key, **** IT IS FIXED HERE ****
//	kp.epos = 21;
	kp.keepShort = false;	// if true, keep short reads; false to discard
	kp.keepN     = false;	// true to keep the reads that have Ns; false to discard

	kp.readx = NULL;
	kp.outprefix = NULL;
	string srstr, nrstr, rpstr;
	int opt;
	//char * optvalue;
	while( (opt = getopt(argc, argv, "i:a:b:o:k:s:n:")) != -1 ) {
		switch( opt ) {
			case 'i': kp.readx=optarg;		break;
			case 'o': kp.outprefix=optarg;	break;
			case 'k': kp.hskip=atoi(optarg);break;
			case 's': srstr=optarg;			break;
			case 'n': nrstr=optarg;			break;
			case '?':
			default : usage( argv[0] );		break;
		}
	}

	// check parameters
	if( kp.readx==NULL || kp.outprefix==NULL ) {
		usage( argv[0] );
	}
	kp.epos = kp.hskip + kp.keylen;
	if( srstr.size() ) {
		if( srstr.compare( "keep" )==0 || srstr.compare( "KEEP" )==0 ) {
			kp.keepShort = true;
		} else if( srstr.compare( "discard" )==0 || srstr.compare( "DISCARD" )==0 ) {
			kp.keepShort = false;
		} else {
			cerr << "Error: invalid -s option, MUST be 'keep' or 'discard'!\n";
			exit(4);
		}
	}
	if( nrstr.size() ) {
		if( nrstr.compare( "keep" )==0 || nrstr.compare( "KEEP" )==0 ) {
			kp.keepN = true;
		} else if( nrstr.compare( "discard" )==0 || nrstr.compare( "DISCARD" )==0 ) {
			kp.keepN = false;
		} else {
			cerr << "Error: invalid -n option, MUST be 'keep' or 'discard'!\n";
			exit(5);
		}
	}
}

int main( int argc, char *argv[] ) {
//	cerr << "Initiate parameters\n";
	kparam kp;
	init_param( kp, argc, argv );

	// prepare files
//	cerr << "Prepare files\n";
	ifstream fqx( kp.readx );
	if( fqx.fail() ) {
		cerr << "Error: read fastq failed!\n";
		fqx.close();
		return 10;
	}

	// container and statistics
	unordered_set<fqKey> log;
	kstat ks;
	ks.uniq = 0;
	ks.dup = 0;
	ks.cN = 0;
	ks.cS = 0;

	// prepare memory
//	cerr << "Prepare memory\n";
	PEREAD *prA, *prB, *wkr, *ldr;	//A, B is real data, wkr and ldr are indicators for work and load
	prA = new PEREAD [ BATCH ];
	prB = new PEREAD [ BATCH ];

//	cerr << "Loading data\n";
	int nwk = load_batch( fqx, prA );
	wkr = prA;
	ldr = prB;

	// start analysis
	bool nextBatch = true;
	int loaded;
	int total = nwk;
	do {
//		cerr << "Working " << total << "\n";
		if( fqx.eof() ) {	//no need to load data, run with 1 thread is enough
			do_rmdup( wkr, nwk, log, ks, kp );
			write_batch( wkr, nwk );
			nextBatch = false;
		} else { // split into 2 threads, 1 for loading data, another for analysis
			omp_set_num_threads( 2 );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				if( tn == 0 ) {	// thread 0: load data
					loaded = load_batch( fqx, ldr );
					total += loaded;
				} else {	// thread 1: analyze data
					do_rmdup( wkr, nwk, log, ks, kp );
					write_batch( wkr, nwk );
				}
			}
			// update parameters for next loop
			nwk = loaded;
			PEREAD *tmp = ldr;
			ldr = wkr;
			wkr = tmp;
		}
	} while( nextBatch );
//	cerr << "Done.\n";
	fqx.close();

	string out = kp.outprefix;
	out += ".log";
	ofstream flog( out.c_str(), ios::app );
	if( flog.fail() ) {
		cerr << "Error: write log failed!\n";
		return 10;
	}
	int valid = ks.uniq + ks.dup;// + cR;
	flog << "Total\t"	<< total
		 << "\nValid\t"	<< valid
		 << "\nUniq\t"	<< ks.uniq
		 << "\nDup\t"	<< ks.dup
		 << "\nShort\t"	<< ks.cS
		 << "\nNs\t"	<< ks.cN
		 << '\n';

	flog.close();

	delete [] prA;
	delete [] prB;
	return 0;
}

