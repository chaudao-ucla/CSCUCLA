/*
 * ALCTChamberPrinter.cpp
 *
 *  Created on: 06 August 2019
 *      Author: Chau Dao
 */

#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <time.h>
#include <stdio.h>
#include <time.h>
#include<map>

#include <TTree.h>
#include <TFile.h>

#include "../include/CSCConstants.h"
#include "../include/CSCClasses.h"
#include "../include/CSCHelperFunctions.h"

#include "../include/CSCInfo.h"
#include "../include/CSCHelper.h"
#include "../include/ALCTHelperFunctions.h"

#include "../include/ALCTChamberPrinter.h"

using namespace std;

int main(int argc, char* argv[]){
	ALCTChamberPrinter p;
	return p.main(argc,argv);
}

int ALCTChamberPrinter::run(string inputfile, unsigned int ST, unsigned int RI, unsigned int CH, unsigned int EC, unsigned int eventnum)
{
	cout << endl << "Running over file: " << inputfile << endl;

	TFile* f = TFile::Open(inputfile.c_str());
	if(!f) throw "Can't open file";

	TTree* t =  (TTree*)f->Get("CSCDigiTree");
	if(!t) throw "Can't find tree";

	ALCTConfig config;

    CSCInfo::Wires wires(t);
	CSCInfo::ALCTs alcts(t);

	/*
	std::vector<ALCTCandidate*> candvec;
	ALCTCandidate * head = new ALCTCandidate(0,1);

	for (int i = 0; i<4; i++)
	{
		ALCTCandidate * cand; 
		cand = (i==0) ? head : new ALCTCandidate(i,1,cand);
	}

	(head->next)->nix();
	head_to_vec(head,candvec);
	for (int i=0; i<candvec.size(); i++)
	{
		cout<<candvec.at(i)<<endl<<endl;
	}*/

	
	for (int e_num = 0; e_num<eventnum; e_num++)
	{
		t->GetEntry(e_num);

		std::vector<ALCT_ChamberHits*> cvec;
		for (int i=0; i<16; i++)
		{
			ALCT_ChamberHits * temp = new ALCT_ChamberHits(ST,RI,CH,EC);
			temp->fill(wires,i);
			cout << *temp << endl; 
			cvec.push_back(temp);
		}
		std::vector<ALCTCandidate*> candvec;

		ALCTCandidate * head = new ALCTCandidate(0,1);

		for (int i = 0; i<(cvec.at(0))->get_maxWi(); i++)
		{
			ALCTCandidate * cand; 
			cand = (i==0) ? head : new ALCTCandidate(i,1,cand);
		}

		preTrigger(cvec,config,head);
		head_to_vec(head,candvec);

		cout << "=== PreTriggerring Results ===" << endl << endl;

		for (int i=0; i<candvec.size(); i++)
		{
			std::cout << candvec.at(i) << endl << endl; 
		}

		patternDetection(cvec, config, head); 
		candvec.clear();
		head_to_vec(head,candvec);

		cout << "=== PatternDetection Results ===" << endl << endl; 

		for (int i=0; i<candvec.size(); i++)
		{
			std::cout << candvec.at(i) << endl << endl; 
		}

		ghostBuster(head);
		clean(head);
		candvec.clear();
		head_to_vec(head,candvec);

		cout << "=== GhostBuster Results ===" << endl << endl; 

		for (int i=0; i<candvec.size(); i++)
		{
			std::cout << candvec.at(i) << endl << endl; 
		}

		cout << " finished emulation results" << endl << endl;

		for (int i=0; i<alcts.size(); i++)
		{
			int chSid1 = CSCHelper::serialize(ST, RI, CH, EC);
			if(!CSCHelper::isValidChamber(ST,RI,CH,EC)) continue;
			int chSid2 = chSid1;
			
			bool me11a	= ST == 1 && RI == 4;
			bool me11b	= ST == 1 && RI == 1;

			if (me11a || me11b)
			{
				if (me11a) chSid2 = CSCHelper::serialize(ST, 1, CH, EC);
				if (me11b) chSid2 = CSCHelper::serialize(ST, 4, CH, EC);
			}
			if (chSid1!=alcts.ch_id->at(i) && chSid2!= alcts.ch_id->at(i)) continue;
			CSCHelper::ChamberId c = CSCHelper::unserialize(alcts.ch_id->at(i));
			cout<< "key wire group = " << (int)(alcts.keyWG->at(i))
				<< ", quality = " << (int)(alcts.quality->at(i)) 
				<< ", accelerator = " << (int) alcts.accelerator->at(i)
				<< ", BX = " <<  (int) alcts.BX->at(i)
				<< ", track number = " << (int) alcts.trkNumber->at(i)
				<< ", Station = " << (int) c.station
				<< ", Ring = " << (int) c.ring
				<< endl; 
		}

		cout << "finished alct for event " << e_num << endl << endl;  

		//if (candvec.size()) return 0; 
	}
    return 0;
}