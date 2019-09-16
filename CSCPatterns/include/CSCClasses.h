/*
 * PatternFitter.h
 *
 *  Created on: Oct 9, 2017
 *      Author: williamnash
 */

#ifndef PATTERNFITTER_H_
#define PATTERNFITTER_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <bitset>
#include <string>
#include <iomanip>

//lambda
#include <functional>

#include "../include/CSCConstants.h"
#include "../include/LUTClasses.h"
#include "../include/CSCInfo.h"

using namespace std;

/* @brief With clever scheme here:
 *
 * https://indico.cern.ch/event/712513/contributions/2959915/attachments/1630667/2599369/20180410_otmb_pattern_logic.pdf
 *
 * Can encode "comparator codes" within patterns.
 *
 */
class ComparatorCode{
public:
	ComparatorCode(bool hits[NLAYERS][3]);
	ComparatorCode(unsigned int comparatorCode);
	ComparatorCode(const ComparatorCode& c);
	ComparatorCode();

	~ComparatorCode() {}

	bool _hits[NLAYERS][3]; //all the hits within a pattern that a chamber scan matched

	unsigned int getLayersMatched() const;
	int getId() const;
	void printCode() const;
	static string getStringInBase4(int code);
	static bool getHits(const unsigned int comparatorCode,bool hits[NLAYERS][3]);

private:
	int _id;
	int _layersMatched;

	void calculateId();
	void calculateLayersMatched();
};




/* @brief CSCPattern is a pattern used by the TMB as a first pass
 * towards creating a CLCT. Each have an acceptance for numerous tracks,
 * but give rough information as to the tracks slope.
 *
 * Comparator codes are created out of hits within a pattern.
 */
class CSCPattern {
public:
	CSCPattern(unsigned int id, bool isLegacy, const bool pat[MAX_PATTERN_WIDTH][NLAYERS]);
	CSCPattern(string name, unsigned int id, bool isLegacy, const bool pat[MAX_PATTERN_WIDTH][NLAYERS]);
	CSCPattern(const CSCPattern &obj);
	CSCPattern();

	~CSCPattern() {}

	const unsigned int _id;
	const bool _isLegacy;
	bool _pat[MAX_PATTERN_WIDTH][NLAYERS];

	CSCPattern makeFlipped(unsigned int id) const;
	CSCPattern makeFlipped(string name, unsigned int id) const;
	void printCode(int code) const;
	int recoverPatternCCCombination(int code, int code_hits[MAX_PATTERN_WIDTH][NLAYERS]) const;
	string getName() const;

private:
	string _name;
};


/* @brief A candidate "road" within the CSC chamber made of a series
 * of hits compressed into a comparator code
 */
class CLCTCandidate {
public:
	CLCTCandidate(CSCPattern p, int horInd, int startTime, bool hits[NLAYERS][3]);
	CLCTCandidate(CSCPattern p, ComparatorCode c, int horInd, int startTime);
	CLCTCandidate(CSCPattern p, int horInd, int startTime, int layMatCount);

	~CLCTCandidate();

	const CSCPattern _pattern;
	const int _horizontalIndex; //half strips, leftmost index of the pattern
	const int _startTime;
	//pointer to whatever LUT Entry is associated with this candidate
	const LUTEntry* _lutEntry;

	int getHits(int code_hits[MAX_PATTERN_WIDTH][NLAYERS]) const;
	float keyStrip() const;
	int keyHalfStrip() const;
	float position() const;
	float slope() const;
	void print3x6Pattern() const;
	void printCodeInPattern() const;
	int comparatorCodeId() const;
	int layerCount() const;
	const ComparatorCode getComparatorCode() const;
	int patternId() const;
	const LUTKey key() const;

	//declare a function pointer used to sort the CLCT Candidates
	typedef function<bool(CLCTCandidate*, CLCTCandidate*)> QUALITY_SORT;
	static QUALITY_SORT quality;

private:
	ComparatorCode* _code;
	int _layerMatchCount;
	//float _quality;

};

class ALCTCandidate
{
	public:
		ALCTCandidate(unsigned int kwg, int pattern);
		ALCTCandidate(unsigned int kwg, int pattern, ALCTCandidate * pred);

		//~ALCTCandidate();

		ALCTCandidate* next;
		ALCTCandidate* prev;

		int get_pattern() const {return _pattern;}
		bool isValid() const {return _isValid;}
		unsigned int get_kwg() const {return _kwg;}
		unsigned int get_first_bx() const {return _first_bx;}
		int get_quality() const {return _quality;}
		int get_tracknumber() const {return _tracknumber;}
		unsigned int get_first_bx_corr() const {return _first_bx_corr;}

		void set_kwg(unsigned int kwg) {_kwg = kwg;}
		void set_first_bx(unsigned int first_bx) {_first_bx = first_bx;}
		void set_first_bx_corr (unsigned int first_bx_corr) {_first_bx_corr = first_bx_corr;}
		void set_quality(unsigned int quality){_quality = quality;}
		void set_pattern(int pattern) {_pattern = pattern;}
		void set_tracknumber(int tracknumber){_tracknumber = tracknumber;}
		
		void flag() {_isValid = false;}
		void nix();

		friend ostream& operator<<(ostream& os, const ALCTCandidate &c);
		friend ostream& operator<<(ostream& os, const ALCTCandidate * const c);

	private:
		int _pattern; 
		int _tracknumber; 
		unsigned int _kwg; 
		unsigned int _first_bx;
		unsigned int _quality;
		bool _isValid;
		unsigned int _first_bx_corr;
};

class ALCTCandidateCollection
{
	public:
		std::vector<int> ch_id; 
		std::vector<int> kwg;
		std::vector<int> first_bx;
		std::vector<int> first_bx_corr;
		std::vector<int> quality; 
		std::vector<int> pattern;
		std::vector<int> track_number;

	ALCTCandidateCollection(TTree *t, int i)
	{
		t->Branch("emulation_ch_id", &ch_id);
		t->Branch("emulation_kwg", &kwg);
		t->Branch("enumation_first_bx", &first_bx);
		t->Branch("emulation_first_bx_corr", &first_bx_corr);
		t->Branch("emulation_quality", &quality);
		t->Branch("emulation_pattern", &pattern);
		t->Branch("emulation_track_number", &track_number);
	}

	void Fill(std::vector<ALCTCandidate*> emulatedALCTs, unsigned int chamberHash);
	void Erase(); 
}

/* @brief Encapsulates hit information for recorded event
 * in a chamber, identified by its station, ring, endcap and chamber
 */
class ChamberHits {
public:
	ChamberHits(unsigned int station=0, unsigned int ring=0,
			unsigned int endcap=0, unsigned int chamber=0, bool isComparator=true);
	ChamberHits(const ChamberHits& c);

	~ChamberHits(){}

	const bool _isComparator;
	const unsigned int _station;
	const unsigned int _ring;
	const unsigned int _endcap;
	const unsigned int _chamber;
	unsigned int minHs() const {return _minHs;}
	unsigned int maxHs() const {return _maxHs;}
	unsigned int nhits() const {return _nhits;}
	float hitMeanHS();
	float hitStdHS();
	int _hits[N_MAX_HALF_STRIPS][NLAYERS];

	bool shift(unsigned int lay) const;

	int fill(const CSCInfo::Comparators& c);
	int fill(const CSCInfo::RecHits& r);
	void print() const; //deprecated
	
	friend ostream& operator<<(ostream& os, const ChamberHits& c);
	ChamberHits& operator-=(const CLCTCandidate& mi);
private:
	unsigned int _nhits;
	unsigned int _minHs; //Halfstrip
	unsigned int _maxHs;

	float _meanHS;
	float _stdHS;
};

class ALCT_ChamberHits
{
	public:
		ALCT_ChamberHits(unsigned int station=0, unsigned int ring=0,
			unsigned int endcap=0, unsigned int chamber=0, bool isWire=true, bool empty = true);
		
		ALCT_ChamberHits(const ALCT_ChamberHits &c);
		
		~ALCT_ChamberHits(){}

		const bool _isWire;
		const unsigned int _station;
		const unsigned int _ring;
		const unsigned int _endcap;
		const unsigned int _chamber;

		unsigned int get_minWi() const {return _minWi;}
		unsigned int get_maxWi() const {return _maxWi;}
		unsigned int get_nhits() const {return _nhits;}
		bool isEmpty() const {return _empty;}
		void regHit() {_empty = false;}

		float get_hitMeanWi();
		float get_hitStdWi();
		
		int _hits[N_KEY_WIRE_GROUPS][NLAYERS];

		ALCT_ChamberHits* prev = 0;
		ALCT_ChamberHits* next = 0; 

		void fill(const CSCInfo::Wires &w);
		void fill(const CSCInfo::Wires &w, int time);
		//void fill(const CSCInfo::Wires &w, int start, int end, int p_ext=6);
		
		friend ostream& operator<<(ostream& os, const ALCT_ChamberHits& c);
		//ALCT_ChamberHits& operator-=(const ALCTCandidate &mi);

	private:
		unsigned int _minWi;
		unsigned int _maxWi;
		unsigned int _nhits;

		bool _empty;

		float _meanWi;
		float _stdWi;
};


#endif /* PATTERNFITTER_H_ */
