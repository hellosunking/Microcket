#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "pairutil.h"
#include "flash2pairs.h"
#include "unc2pairs.h"

using namespace std;

/*
 * Author: Kun Sun @ SZBL
*/

int main( int argc, char *argv[] ) {
	if( argc < 4 ) {
		cerr << "\nUsage: " << argv[0] << " <in.sam> <mode=flash|unc> <out.prefix> [thread=4]"
			 << "\n\nTask: extract the pairs from the alignment result."
			 << "\n3 files will be written: out.mode.pairs out.mode.stat out.mode.sam"
			 << "\n\nThis program is part of Microcket, and is NOT supposed to be called manually by the user.\n\n";
		exit(2);
	}

	int thread = 4;
	if( argc > 4 ) {
		thread = atoi( argv[4] );
		if( thread < 2 ) {
			cerr << "Error: at least 2 threads are required.\n";
			return 5;
		}
	}
//	cerr << "INFO: " << thread << " threads will be used.\n";
	void (*sam2pair)(vector<string> *, vector<string> *, string &, string &, kstat &);

	string base = argv[2];
	if( base.compare("flash") == 0 ) {
		sam2pair = flash2pairs;
	} else if( base.compare("unc") == 0 ) {
		sam2pair = unc2pairs;
	} else {
		cerr << "Error: Unknown mode, must be 'flash' or 'unc'.\n";
		return 6;
	}

//	cerr << "Prepare files ...\n";
	ifstream fin( argv[1] );
	if( fin.fail() ) {
		cerr << "Error: read input file failed!\n";
		fin.close();
		return 10;
	}

	base = argv[3];
	base += '.';
	base += argv[2];
	string out = base;
	out += ".sam";
	FILE *fsam = fopen( out.c_str(), "wb");
	if( fsam==NULL ) {
		cerr << "Error: write output file failed!\n";
		fin.close();
		fclose( fsam );
		return 11;
	}

//	cerr << "Prepare stat and memeory\n";
	// container and statistics
	kstat *ks = new kstat [ thread ];
	for( int i=0; i!=thread; ++i ) {
		ks[i].lowMap    = 0;
		ks[i].manyHits  = 0;
		ks[i].unpaired  = 0;
		ks[i].selfCircle= 0;
		ks[i].trans     = 0;
		ks[i].cis0      = 0;
		ks[i].cis1K     = 0;
		ks[i].cis10K    = 0;
	}

	// prepare memory
	vector<string> *prA = new vector<string>[ BATCH ];
	vector<string> *prB = new vector<string>[ BATCH ];
	string line;
	// load the first line
	stringstream ss;
	string unk;
	unsigned int mapinfo, mapQ;
//	cerr << "Load the 1st batch\n";
	while( getline(fin, line) ) {
		if( line[0] == '@') {	// sam header
			continue;
		}

		ss.str( line );
		ss >> unk >> mapinfo >> unk >> unk >> mapQ;
		// discard 2nd alignment, QC failed, PCR duplicate, or low mapQ
		if( mapinfo & 0x700 || mapQ < min_mapQ ) {
			continue;
		}

		// meet a valid record
		break;
	}
	int loaded = load_batch( fin, prA, line, BATCH );
//	cerr << loaded << " records loaded\n";

	// result containers
	string *outPair = new string [ BATCH ];
	string *outSam  = new string [ BATCH ];

	// start analysis
	vector<string> *wkr = prA;
	vector<string> *ldr = prB;
	bool nextBatch = true;
	do {
		if( fin.eof() ) {	//no need to load data, use all threads to do the analysis
			omp_set_num_threads( thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				unsigned int start = loaded *  tn    / thread;
				unsigned int end   = loaded * (tn+1) / thread;
				(*sam2pair)( wkr+start, wkr+end, outPair[tn], outSam[tn], ks[tn] );

				fwrite( outPair[tn].c_str(), 1, outPair[tn].size(), stdout );
				fwrite(  outSam[tn].c_str(), 1, outSam[tn].size(),  fsam  );
				outPair[tn].clear();
				outSam[tn].clear();
			}
			nextBatch = false;
		} else { // split into N threads, 1 for loading data, others for analysis
			omp_set_num_threads( thread );
			#pragma omp parallel
			{
				unsigned int tn = omp_get_thread_num();
				unsigned int wkt = thread - 1;
				if( tn == wkt ) {	//last thread: load data
					loaded = load_batch( fin, ldr, line, BATCH );
				} else {	// others: analyze data
					unsigned int start = loaded *  tn    / wkt;
					unsigned int end   = loaded * (tn+1) / wkt;
					(*sam2pair)( wkr+start, wkr+end, outPair[tn], outSam[tn], ks[tn] );
					fwrite( outPair[tn].c_str(), 1, outPair[tn].size(), stdout );
					fwrite(  outSam[tn].c_str(), 1, outSam[tn].size(),  fsam  );
					outPair[tn].clear();
					outSam[tn].clear();
				}
			}
			// update parameters for next loop
			vector<string> *tmp = ldr;
			ldr = wkr;
			wkr = tmp;
		}
//		cerr << "Write output\n";
	} while( nextBatch );
//	cerr << "Done.\n";
	fin.close();
	fclose( fsam );

	out = base;
	out += "2pairs.log";
	ofstream fst( out.c_str() );
	if( fst.fail() ) {
		cerr << "Error: write log file failed!\n";
		return 10;
	}
	for( int i=1; i!=thread; ++i ) {
		ks[0].lowMap  += ks[i].lowMap;
		ks[0].manyHits+= ks[i].manyHits;
		ks[0].unpaired+= ks[i].unpaired;
		ks[0].trans += ks[i].trans;
		ks[0].cis0  += ks[i].cis0;
		ks[0].cis1K += ks[i].cis1K;
		ks[0].cis10K+= ks[i].cis10K;
	}
	fst << "lowMap\t" << ks[0].lowMap
		<< "\nmanyHits\t" << ks[0].manyHits
		<< "\nunpaired\t" << ks[0].unpaired
		<< "\nselfCircle\t" << ks[0].selfCircle
		<< "\ntrans\t" << ks[0].trans
		<< "\ncis10K\t" << ks[0].cis10K
		<< "\ncis1K\t" << ks[0].cis1K
		<< "\ncis0\t" << ks[0].cis0 << '\n';
	fst.close();

	delete [] prA;
	delete [] prB;
	delete [] ks;
	delete [] outPair;
	delete [] outSam;

	return 0;
}

