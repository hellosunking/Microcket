#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "pairutil.h"

using namespace std;

/*
 * Author: Kun Sun @ SZBL
*/

extern bool writeSam;

void flash2pairs( vector<string> *start, vector<string> *end, string &outPair, string &outSam, kstat &ks ) {
	stringstream ss;
	string rid, chr1, chr2, cigar1, cigar2, mapQ;
	unsigned int mapinfo, pos1, pos2;
	char strand1, strand2;
	SEGMENT s1, s2;
	string p1s, p2s;

	for( vector<string> *pvs=start; pvs!=end; ++pvs ) {
		if( pvs->size() == 1 ) {	// only 1 record, could be xxxM, or xxMyyNzzM
			ss.str( pvs->at(0) );
			ss.clear();
			ss >> rid >> mapinfo >> chr1 >> pos1 >> mapQ >> cigar1;
//			cerr << rid << " 1 record\n";
//			s1.strand = (mapinfo & 16) ? '-' : '+';
			cigar2segment(cigar1, pos1, s1);

			if( s1.segCnt > 2 ) {	// tow many segments. TODO: rescue 100M200N30M10000N100M thing?
//				cerr << "E2\t" << rid << " has too many introns.\n";
				++ ks.manyHits;
				continue;
			}
			// 1-2 segment, chould be xxSyyMzzS, xxMyyNzzM
			//check read integrity
			if( ! check_integrity_1_seg(s1) ) {
//				cerr << "E1\t" << rid << " low map.\n";
				++ ks.lowMap;
				continue;
			}

			// output the outmost loci
			pos2 = s1.right[ s1.segCnt-1 ];
			unsigned int dist = pos2 - pos1;
//			if( dist <= maxSelfCircleDist ) {	// this should not be possible?
//				++ ks.selfCircle;	// self-loop
//				continue;
//			}
			if(dist>=10000) {++ ks.cis10K;} else if (dist>=1000) {++ ks.cis1K;} else {++ ks.cis0;}

			ss.clear();
			ss.str("");
			ss << pos1 << '\t' << pos2;
			ss >> p1s >> p2s;
			outPair += rid;outPair+='\t';
			outPair += chr1;outPair+='\t';outPair += p1s;outPair+='\t';
			outPair += chr1;outPair+='\t';outPair += p2s;
			outPair += "\t+\t-\n";

			if( writeSam ) {
				outSam += pvs->at(0);
				outSam += '\n';
			}
		} else if( pvs->size() == 2 ) {	// two segments, common scenario
			ss.str( pvs->at(0) );
			ss.clear();
			ss >> rid >> mapinfo >> chr1 >> pos1 >> mapQ >> cigar1;
			strand1 = (mapinfo & 16) ? '-' : '+';
			cigar2segment(cigar1, pos1, s1);

			ss.str( pvs->at(1) );
			ss.clear();
			ss >> rid >> mapinfo >> chr2 >> pos2 >> mapQ >> cigar2;
			strand2 = (mapinfo & 16) ? '-' : '+';
			cigar2segment(cigar2, pos2, s2);

			if( s1.segCnt!=1 || s2.segCnt!=1 ) {
//				cerr << "E2\t" << rid << " has 2 hits AND introns!\n";
				++ ks.manyHits;
				continue;
			}

			// check read integrity
//			cerr << "Check integrity\n";
			if( ! check_integrity_2_seg(s1, s2) ) {
//				cerr << "E1\t" << rid << " low map.\n";
				++ ks.lowMap;
				continue;
			}

			//get the outmost part
			if( s1.leftClip > s1.rightClip ) {	//segment is clipped on the left, use the right part
				pos1 = s1.right[0];
			}
			if( s2.leftClip > s2.rightClip ) {	//segment is clipped on the left, use the right part
				pos2 = s2.right[0];
			}
		
			//output
			ss.clear();
			ss.str("");
			ss << pos1 << '\t' << pos2;
			ss >> p1s >> p2s;

			int chrcmp = chr1.compare(chr2);
			if( (chrcmp<0) || (chrcmp==0 && pos1 < pos2) ) {
				if( chrcmp == 0 ) {
					unsigned int dist = pos2 - pos1;
					if( dist <= maxSelfCircleDist ) {
						++ ks.selfCircle;	// self-loop
						continue;
					}
					if(dist>=10000) {++ ks.cis10K;} else if (dist>=1000) {++ ks.cis1K;} else {++ ks.cis0;}
				} else {
					++ ks.trans;
				}

				outPair += rid; outPair+='\t';
				outPair += chr1;outPair+='\t';outPair += p1s; outPair+='\t';
				outPair += chr2;outPair+='\t';outPair += p2s; outPair+='\t';
				outPair += strand1;outPair+='\t';outPair += strand2;
				outPair += '\n';
			} else {
				if( chrcmp == 0 ) {
					unsigned int dist = pos1 - pos2;
					if( dist <= maxSelfCircleDist ) {
						++ ks.selfCircle;	// self-loop
						continue;
					}
					if(dist>=10000) {++ ks.cis10K;} else if (dist>=1000) {++ ks.cis1K;} else {++ ks.cis0;}
				} else {
					++ ks.trans;
				}
				outPair += rid; outPair+='\t';
				outPair += chr2;outPair+='\t';outPair += p2s; outPair+='\t';
				outPair += chr1;outPair+='\t';outPair += p1s; outPair+='\t';
				outPair += strand2;outPair+='\t';outPair += strand1;
				outPair += '\n';
			}

			if( writeSam ) {
				outSam += pvs->at(0);outSam += '\n';
				outSam += pvs->at(1);outSam += '\n';
			}
		} else {
//			cerr << "E9\t" << rid << " has too many segments.\n";
			++ ks.manyHits;
		}
	}
}

