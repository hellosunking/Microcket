#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/*
 * Author: Ahfyth
 */

int main( int argc, char *argv[] ) {
	ios::sync_with_stdio( false );
	if( argc < 4 ) {
		cerr << "\nUsage: " << argv[0] << " <read1.fq> <read2.fq> <output.prefix> [cut=10] [min.size=36]\n\n";
		return 2;
	}
	unsigned int cutSize  = 10;
	unsigned int keepSize = 36;

	if( argc > 4 ) {
		cutSize = atoi( argv[4] );
		if( argc > 5 ) {
			keepSize = atoi( argv[5] );
		}
	}

	// prepare files
	ifstream fq1( argv[1] ), fq2( argv[2] );
	if( fq1.fail() || fq2.fail() ) {
		cerr << "Error: read fastq failed!\n";
		fq1.close();
		fq2.close();
		return 10;
	}
	string base = argv[3];
	string out = base;
	out += ".read1.fq";
	ofstream fout1( out.c_str() );
	out = base;
	out += ".read2.fq";
	ofstream fout2( out.c_str() );
	if( fout1.fail() || fout2.fail() ) {
		cerr << "Error: write output file failed!\n";
		fq1.close();
		fq2.close();
		fout1.close();
		fout2.close();
		return 11;
	}

	// start
	register unsigned int minSize = cutSize + keepSize;
	register unsigned int fail=0, pass=0;

	string id1, seq1, qual1, id2, seq2, qual2, unk;
	while( 1 ) {
		getline( fq1, id1 );
		if( fq1.eof() )break;
		getline( fq1, seq1 );
		getline( fq1, unk );
		getline( fq1, qual1 );

		getline( fq2, id2 );
		getline( fq2, seq2 );
		getline( fq2, unk );
		getline( fq2, qual2 );

		// assume that seq1 is the same size as seq2
		if( seq1.size() < minSize ) {	// too short
			++ fail;
			continue;
		}

		int s = seq1.size() - cutSize;
		seq1.resize( s );
		qual1.resize( s );
		seq2.resize( s );
		qual2.resize( s );

		++ pass;
		fout1 << id1 << '\n' << seq1 << "\n+\n" << qual1 << '\n';
		fout2 << id2 << '\n' << seq2 << "\n+\n" << qual2 << '\n';
	}
	fq1.close();
	fq2.close();
	fout1.close();
	fout2.close();

	out = base;
	out += ".log";
	ofstream fout( out.c_str() );
	if( fout.fail() ) {
		cerr << "Error: write log failed!\n";
		return 10;
	}

	fout<< "Read1\t" << argv[1]
		<< "\nRead2\t" << argv[2]
		<< "\nTotal\t" << pass+fail
		<< "\nPass\t" << pass << '\n';
	fout.close();

	return 0;
}

