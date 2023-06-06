#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
//#include <time.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#ifndef __PAIRUTIL__
#define __PAIRUTIL__

//structures
// statistics
typedef struct {
	unsigned int lowMap;
	unsigned int manyHits;
	unsigned int unpaired;
	unsigned int trans;
	unsigned int selfCircle;
	unsigned int cis0;
	unsigned int cis1K;
	unsigned int cis10K;
} kstat;

// segments
typedef struct {
//	string chr;
//	char strand;
	int segCnt;
	int leftClip;
	int rightClip;
	vector<int> left;
	vector<int> right;
	int mappable;
} SEGMENT;

typedef struct {
	string chr;
	int pos;
	char strand;
	string cigar;
} SEGKEY;

//parameters
const int BATCH = 1 << 18;
const int maxSegNum = 16;
static int min_mapQ = 10;				//minimum mapQ; NOTE: it is changed to a parameter in v1.2
static float min_mapped_ratio = 0.5;	//the read must be covered at least 80% to consider keep; NOTE: it is changed to a parameter in v1.2
static bool writeSam = true;
const int min_clip_size = 20;			//minimum bp to be considered as a clip; STAR sometimes give XXXM1S thing
const int max_unmapped_allowed = 20;	//minimum unmapped bp to discard the read, this is to rescue 100M15S thing

const int maxSelfCircleDist = 100;		//NOTE: changed to 100 in v1.2
const int maxPairDist = 1000;			//for check_pair

//functions
//return the status of this function
//mappable and clips are saved in SEGRESULT
bool cigar2segment(const string &cigar, const int start, SEGMENT &sr) {
//	cerr << "==== Cigar " << cigar;
	// init/reset sr
	sr.segCnt = 0;
	sr.leftClip = 0;
	sr.rightClip= 0;
	sr.left.clear();
	sr.right.clear();
	sr.mappable = 0;

	// init SEGMETS
	int index = 0;
	sr.left.push_back( start );
	sr.right.push_back( 0 );

	int currPos = start;
	int cigarLen= cigar.length();
	int value = 0;
	for( int j=0; j!=cigarLen; ++j ) {
		char c = cigar[j];
		if( c>='0' && c<='9' ) {	// number
			value *= 10;
			value += c-'0';
		} else {
			if( c=='H' || c=='S' ) {	// hard/soft clip, skip and mark unmapped
				if( j==cigarLen-1 ) {
					sr.rightClip = value;
				} else if( index==0 ) {
					sr.leftClip = value;
				} else {
					cerr << "ERROR CLIP!\n";
					return false;
				}
			} else if( c == 'M' ) {
				sr.mappable += value;
				currPos += value;
				sr.right[index] = currPos - 1;
			} else if( c == 'D' ) {	// deletion, move the cursor
				currPos += value;
				sr.right[index] = currPos - 1;	//this is unnecessary as D must be in the middle and followed by xxM
			} else if( c == 'I' ) {	// insertion, skip it
				// do nothing
			} else if( c == 'N' ) {	// N for introns, initiate a new segment
				currPos += value;
				++ index;
				sr.left.push_back( currPos );
				sr.right.push_back( 0 );
			} else {
				cerr << "ERROR: unknown element in " << cigar << "!\n";
				return false;
			}
			value = 0;
		}
	}

	// double-check
	if( sr.right[index] == 0 ) {
		cerr << "ERROR: something is wrong for cigar " << cigar << "!\n";
		return false;
	}
	sr.segCnt = index + 1;
//	cerr << ", " << sr.segCnt << " segments extracted.\n";
	return true;
}

// check if 2 segments could form a valid pair
// NOTE: I only check positions here! chr and strand should be checked outside
//inline bool canPair( const SEGMENT &s1, unsigned int index1, const SEGMENT &s2, unsigned int index2 ) {
//	return s1.left[index1]<s2.left[index2] && s2.right[index2]-s1.left[index1]<=maxPairDist;
//}

//return the number of records loaded
//the line contains the 1st record, and after the function called, it contains the last record (if any)
int load_batch( ifstream &fin, vector<string> *records, string &line, int maxNum ) {
	// clean up
	for( int i=0; i!=maxNum; ++i ) {
		records[i].clear();
	}

	stringstream ss;
	string lastId, currId, mapflag, unk;
	unsigned int mapinfo, mapQ;
	ss.str( line );
	ss.clear();
	ss >> lastId;

	records[0].push_back( line );

	int cnt = 0;
	while( getline(fin, line) ) {
		ss.str( line );
		ss.clear();
		ss >> currId >> mapinfo >> unk >> unk >> mapQ;

		if( mapQ < min_mapQ )	// low quality alignment
			continue;
		
		if( mapinfo & 0x700 )		//0x700=1024+512+256; 256:2nd alignment, 512: QC failed, 1024: PCR duplicate
			continue;

		if( currId != lastId ) {	// a new fragment comes
			++ cnt;
			if( cnt == BATCH ) {
				break;
			} else {
				records[ cnt ].push_back( line );
				lastId = currId;
			}
		} else {	// same fragment as the previous one
			records[ cnt ].push_back( line );
		}
	}

	return cnt;
}

// check read integrity
inline bool check_integrity_1_seg( const SEGMENT &s ) {
	int total = s.mappable;
	if( s.leftClip > min_clip_size )	// this is to save the 5S100M thing
		total += s.leftClip;
	if( s.rightClip > min_clip_size )
		total += s.rightClip;

	return s.mappable >= total*min_mapped_ratio;
}

inline bool check_integrity_2_seg( const SEGMENT &s1, const SEGMENT &s2 ) {
	int total_1 = s1.mappable;
	if( s1.leftClip > min_clip_size )	// this is to save the 5S100M thing
		total_1 += s1.leftClip;
	if( s1.rightClip > min_clip_size )
		total_1 += s1.rightClip;

	int total_2 = s2.mappable;
	if( s2.leftClip > min_clip_size )	// this is to save the 5S100M thing
		total_2 += s2.leftClip;
	if( s1.rightClip > min_clip_size )
		total_2 += s2.rightClip;

	if( total_1 > total_2 ) {
		return s1.mappable + s2.mappable >= total_1*min_mapped_ratio;
	} else {
		return s1.mappable + s2.mappable >= total_2*min_mapped_ratio;
	}
}

#endif

