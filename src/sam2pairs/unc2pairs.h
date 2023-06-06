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

void unc2pairs( vector<string> *start, vector<string> *end, string &outPair, string &outSam, kstat &ks ) {
	stringstream ss;
	string rid, mapQ;
	unsigned int mapinfo, pos1, pos2;
	SEGMENT s1, s2, s3;
	string *chr1, *chr2;
	char strand1, strand2;
	string p1s, p2s;
	SEGKEY segkey;
	vector<SEGKEY> R1, R2;

	int category;	//0: 1+1, 1: 1+2, 2: 2+1

	for( vector<string> *pvs=start; pvs!=end; ++pvs ) {
		// separate Read1 and Read2
		R1.clear();
		R2.clear();
		for( int i=0; i!=pvs->size(); ++i ) {
		    ss.str( pvs->at(i) );
			ss.clear();
			ss >> rid >> mapinfo >> segkey.chr >> segkey.pos >> mapQ >> segkey.cigar;
//			cerr << rid << "\n";
			segkey.strand = (mapinfo & 16 ) ? '-' : '+';
			if( mapinfo & 64 ) {
				R1.push_back( segkey );
//				cerr << "R1\t" << segkey.cigar << "\n";
			} else if( mapinfo & 128 ) {
				R2.push_back( segkey );
//				cerr << "R2\t" << segkey.cigar << "\n";
			} else {	//unknown flag, discard it
//				cerr << "E0\t" << rid << " unacceptable mapinfo " << mapinfo << "!\n";
				continue;
			}
		}

		// check segment numbers, currently support 1 + 1/2 combination
		if( R1.size()==0 || R2.size()==0 ) {	//something went wrong, discard; TODO: maybe I can rescue one?
//			cerr << "E1\t" << rid << " unpaired.\n";
			continue;
		}
		if( R1.size() + R2.size() > 3 ) {	// either read has too many hits
//			cerr << "E1\t" << rid << " too many hits.\n";
			continue;
		}

		if( R1.size() == 1 ) {	// only 1 record, could be xxxM, or xxMyyNzzM
			cigar2segment(R1.at(0).cigar, R1.at(0).pos, s1);
			//check read integrity
			if( ! check_integrity_1_seg(s1) ) {
				++ ks.lowMap;
				continue;
			}

			if( R2.size() == 1 ) {
				cigar2segment(R2.at(0).cigar, R2.at(0).pos, s2);
				//check read integrity
				if( ! check_integrity_1_seg(s2) ) {
					++ ks.lowMap;
					continue;
				}

				if( s1.segCnt + s2.segCnt > 3 ) {	// tow many segments. TODO: rescue 100M200N30M10000N100M thing?
//					cerr << "E2\t" << rid << " has too many introns.\n";
					++ ks.manyHits;
					continue;
				}

				category = 0;
			} else {
				cigar2segment(R2.at(0).cigar, R2.at(0).pos, s2);
				cigar2segment(R2.at(1).cigar, R2.at(1).pos, s3);
				//check read integrity
				if( ! check_integrity_2_seg(s2, s3) ) {
					++ ks.lowMap;
					continue;
				}
				if( s1.segCnt!=1 || s2.segCnt!=1 || s3.segCnt!=1 ) {
//					cerr << "E1\t" << rid << " too many hits.\n";
					++ ks.manyHits;
					continue;
				}
				category = 1;
			}
		} else {
			cigar2segment(R1.at(0).cigar, R1.at(0).pos, s1);
			cigar2segment(R1.at(1).cigar, R1.at(1).pos, s2);
			//check read integrity
			if( ! check_integrity_2_seg(s1, s2) ) {
				++ ks.lowMap;
				continue;
			}

			cigar2segment(R2.at(0).cigar, R2.at(0).pos, s3);
			//check read integrity
			if( ! check_integrity_1_seg(s3) ) {
				++ ks.lowMap;
				continue;
			}

			if( s1.segCnt!=1 || s2.segCnt!=1 || s3.segCnt!=1 ) {
//				cerr << "E1\t" << rid << " too many hits.\n";
				++ ks.manyHits;
				continue;
			}
			category = 2;
		}

		// deal with each category of segments
		// get the chr1/pos1/strand1, chr2/pos2/strand2
		if( category == 0 ) {	// 1+1, need to check 100M200N30M10000N100M thing
			strand1 = R1.at(0).strand;
			strand2 = R2.at(0).strand;
			chr1 = & R1.at(0).chr;
			chr2 = & R2.at(0).chr;

			// 1-2 segment, chould be xxSyyMzzS, xxMyyNzzM
			// NOTE: Read1 and Read2 MUST have pairable segments if any of them have > 1 segments
			if( s1.segCnt==1 && s2.segCnt==1 ) {	//both reads only have only 1 segment, simple scenario
				// get the outmost of seg1, whose cigar could be 100M, or 5S95M
				if( strand1 == '+' ) {
					pos1 = s1.left[0]; //use leftmost
				} else {
					pos1 = s1.right[0];	//use rightmost
				}

				if( strand2 == '+' ) {
					pos2 = s2.left[0]; //use leftmost
				} else {
					pos2 = s2.right[0];	//use rightmost
				}
			} else if( s1.segCnt==2 ) {	// s2.segCnt must be 1; check whether they are pairable
				if( strand1 == '+' ) {	// R1 ====-------=====.....======== R2
					if( strand2=='-' && chr1->compare(*chr2)==0 &&
							s1.left[1]<s2.left[0] && s2.right[0]-s1.left[1]<=maxPairDist ) {
						pos1 = s1.left[0];
						pos2 = s2.right[0];
					} else {
//						cerr << "E10\t" << rid << " unpairable!\n";
						++ ks.unpaired;
						continue;
					}
				} else {	// R2 ========.....=====--------=== R1
					if( strand2=='+' && chr1->compare(*chr2)==0 &&
							s2.left[0]<s1.left[0] && s1.right[0]-s2.left[0]<=maxPairDist ) {
						pos1 = s1.right[1];
						pos2 = s2.left[0];
					} else {
//						cerr << "E10\t" << rid << " unpairable!\n";
						++ ks.unpaired;
						continue;
					}
				}
			} else {	//s1.segCnt==1 and s2.segCnt==2
				if( strand1 == '+' ) {	// R1 ========...=====-----=== R2
					if( strand2=='-' && chr1->compare(*chr2)==0 &&
							s1.left[0]<s2.left[0] && s2.right[0]-s1.left[0]<=maxPairDist ) {
						pos1 = s1.left[0];
						pos2 = s2.right[1];
					} else {
//						cerr << "E10\t" << rid << " unpairable!\n";
						++ ks.unpaired;
						continue;
					}
				} else {	// R2 ===----=====....======== R1
					if( strand2=='+' && chr1->compare(*chr2)==0 &&
							s2.left[1]<s1.left[0] && s1.right[0]-s2.left[1]<=maxPairDist ) {
						pos1 = s1.right[0];
						pos2 = s2.left[0];
					} else {
//						cerr << "E10\t" << rid << " unpairable!\n";
						++ ks.unpaired;
						continue;
					}
				}
			}
		} else if( category == 1 ) {	// 1+2
			//check whether R1 (s1) could pair to R2's 1 segment (s2 or s3)
			strand1 = R1.at(0).strand;
			chr1 = & R1.at(0).chr;

			int mate = 0;
			// check s1 vs s2
			if( strand1 == '+' ) {	// s1 ========...=== s2 / ======= s3
				if( R2.at(0).strand=='-' && chr1->compare(R2.at(0).chr)==0 &&
						s1.left[0]<s2.left[0] && s2.right[0]-s1.left[0]<=maxPairDist ) {
					pos1 = s1.left[0];
					mate = 2;
					// chr2/pos2/strand2 should be extracted from s3
				}
			} else {	// s2 ========....======== s1
				if( R2.at(0).strand=='+' && chr1->compare(R2.at(0).chr)==0 &&
							s2.left[0]<s1.left[0] && s1.right[0]-s2.left[0]<=maxPairDist ) {
					pos1 = s1.right[0];
					mate = 2;
				}
			}

			if( mate == 0 ) {	// s1 and s2 could not pair, check s1 vs s3
				if( strand1 == '+' ) {	// s1 ========...=== s3 / ======= s2
					if( R2.at(1).strand=='-' && chr1->compare(R2.at(1).chr)==0 &&
							s1.left[0]<s3.left[0] && s3.right[0]-s1.left[0]<=maxPairDist ) {
						pos1 = s1.left[0];
						mate = 3;
					}
				} else {	// s3 ========....======== s1
					if( R2.at(1).strand=='+' && chr1->compare(R2.at(1).chr)==0 &&
								s3.left[0]<s1.left[0] && s1.right[0]-s3.left[0]<=maxPairDist ) {
						pos1 = s1.right[0];
						mate = 3;
					}
				}
			}

			if( mate == 0 ) {	// unpairable
//				cerr << "E10\t" << rid << " unpairable!\n";
				++ ks.unpaired;
				continue;
			} else if( mate == 2 ) {	// s1 and s2 are pairable, chr2/pos2/strand2 should be extracted from s3
				chr2 = & R2.at(1).chr;
				strand2 = R2.at(1).strand;
				if( s3.leftClip > s3.rightClip ) {	// clip is on the left, use the rightmost
					pos2 = s3.right[0];
				} else {
					pos2 = s3.left[0];
				}
			} else {
				chr2 = & R2.at(0).chr;
				strand2 = R2.at(0).strand;
				if( s2.leftClip > s2.rightClip ) {	// clip is on the left, use the rightmost
					pos2 = s2.right[0];
				} else {
					pos2 = s2.left[0];
				}
			}
		} else {	// 2+1
			//check whether R1's 1 segment (s1 or s2) could pair to R2 (s3)
			strand2 = R2.at(0).strand;
			chr2 = & R2.at(0).chr;

			int mate = 0;
			// check s1 vs s3
			if( strand2 == '+' ) {	// s3 ========....======= s1
				if( R1.at(0).strand=='-' && chr2->compare(R1.at(0).chr)==0 &&
						s3.left[0]<s1.left[0] && s1.right[0]-s3.left[0]<=maxPairDist ) {
					pos2 = s3.left[0];
					mate = 1;
				}
			} else {	// s1 ========....======== s3
				if( R1.at(0).strand=='+' && chr2->compare(R1.at(0).chr)==0 &&
						s1.left[0]<s3.left[0] && s3.right[0]-s1.left[0]<=maxPairDist ) {
					pos2 = s3.right[0];
					mate = 1;
				}
			}

			if( mate == 0 ) {	// s1 and s3 could not pair, check s2 vs s3
				if( strand2 == '+' ) {	// s3 ========....======= s2
					if( R1.at(1).strand=='-' && chr2->compare(R1.at(1).chr)==0 &&
							s3.left[0]<s2.left[0] && s2.right[0]-s3.left[0]<=maxPairDist ) {
						pos2 = s3.left[0];
						mate = 2;
					}
				} else {	// s2 ========....======== s3
					if( R1.at(1).strand=='+' && chr2->compare(R1.at(1).chr)==0 &&
							s2.left[0]<s3.left[0] && s3.right[0]-s2.left[0]<=maxPairDist ) {
						pos2 = s3.right[0];
						mate = 2;
					}
				}
			}

			if( mate == 0 ) {	// unpairable
//				cerr << "E10\t" << rid << " unpairable!\n";
				++ ks.unpaired;
				continue;
			} else if( mate == 1 ) {	// s1 and s3 are pairable, chr1/pos1/strand1 should be extracted from s2
				chr1 = & R1.at(1).chr;
				strand1 = R1.at(1).strand;
				if( s2.leftClip > s2.rightClip ) {	// clip is on the left, use the rightmost
					pos1 = s2.right[0];
				} else {
					pos1 = s2.left[0];
				}
			} else {
				chr1 = & R1.at(0).chr;
				strand1 = R1.at(0).strand;
				if( s1.leftClip > s1.rightClip ) {	// clip is on the left, use the rightmost
					pos1 = s1.right[0];
				} else {
					pos1 = s1.left[0];
				}
			}
		}

		// output
		ss.clear();
		ss.str("");
		ss << pos1 << '\t' << pos2;
		ss >> p1s >> p2s;
		int chrcmp = chr1->compare( *chr2 );
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
			outPair += *chr1;outPair+='\t'; outPair += p1s;outPair+='\t';
			outPair += *chr2;outPair+='\t'; outPair += p2s;outPair+='\t';
			outPair += strand1;outPair+='\t'; outPair += strand2;
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
			outPair += *chr2;outPair+='\t'; outPair += p2s;outPair+='\t';
			outPair += *chr1;outPair+='\t'; outPair += p1s;outPair+='\t';
			outPair += strand2;outPair+='\t'; outPair += strand1;
			outPair += '\n';
		}

		// deal with sam
		if( writeSam ) {
			for( int i=0; i!=pvs->size(); ++i ) {
				outSam += pvs->at(i);
				outSam += '\n';
			}
		}
	}
}

