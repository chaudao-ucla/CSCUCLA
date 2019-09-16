/*
 * ALCTEmulationTreeCreator.cpp
 *
 *  Created on: 29 July 2019
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
#include <map>

#include <TTree.h>
#include <TFile.h>

#include "../include/CSCConstants.h"
#include "../include/CSCClasses.h"
#include "../include/CSCHelperFunctions.h"
#include "../include/ALCTHelperFunctions.h"

#include "../include/CSCInfo.h"
#include "../include/CSCHelper.h"

#include "../include/ALCTEmulationTreeCreator.h"

using namespace std; 

int main(int argc, char* argv[])
{
    ALCTEmulationTreeCreator p;
    return p.main(argc, argv);
}

int ALCTEmulationTreeCreator::run(string inputfile, string outputfile, int start, int end)
{
	/**********************
	 * TREE INITIALIZATION
	 **********************/

    auto t1 = std::chrono::high_resolution_clock::now();

    cout << endl << "Running over file: " << inputfile << endl;

    TFile* f = TFile::Open(inputfile.c_str());
    if (!f) throw "Can't open file";

    TTree* t = (TTree*) f->Get("CSCDigiTree");
    if(!t) throw "Can't find Tree";

    TFile* outF = new TFile(outputfile.c_str(),"RECREATE");
	if(!outF)
	{
		printf("Failed to open output file: %s\n", outputfile.c_str());
		return -1;
	}

    TTree* alct_t_emu = new TTree("ALCTEmulationResults", "ALCTEmulationResults");

    CSCInfo::ALCTs alcts(t);
	CSCInfo::LCTs lcts(t);
    CSCInfo::Wires wires(t);
	CSCInfo::Segments segments(t); 

    /**********************
	 * EVENT LOOP
	 **********************/

	if(end > t->GetEntries() || end < 0) end = t->GetEntries();

	cout << endl;
	
	cout << "Starting Event: " << start << " Ending Event: " << end << endl << endl;

	ALCTConfig config;
	//config.set_narrow_mask_flag(true);

	for(int i = start; i < end; i++) 
	{
		if(!(i%100)) printf("%3.2f%% Done --- Processed %u Events\n\n", 100.*(i-start)/(end-start), i-start);

		t->GetEntry(i);

		/**********************
	 	* CHAMBER LOOP
	 	**********************/

		for(unsigned int chamberHash = 0; chamberHash < (unsigned int)CSCHelper::MAX_CHAMBER_HASH; chamberHash++)
		{
			CSCHelper::ChamberId c = CSCHelper::unserialize(chamberHash);

			unsigned int EC = c.endcap;
			unsigned int ST = c.station;
			unsigned int RI = c.ring;
			unsigned int CH = c.chamber;


			if (ST == 1 && RI == 1) continue; 	// we skip these chambers because the degeneracy is accounted for in fill functions
												// in ALCTChambers::fill()

			if(!CSCHelper::isValidChamber(ST,RI,CH,EC)) continue;

			/**********************
	 		* CHAMBER LIST SETUP
	 		**********************/
			
			std::vector<ALCT_ChamberHits*> cvec;

			for (int i=0; i<16; i++)
			{
				ALCT_ChamberHits * temp = new ALCT_ChamberHits(ST,RI,CH,EC);
				temp->fill(wires,i);
				cvec.push_back(temp);
			}

			std::vector<std::vector<ALCTCandidate*>> end_vec; 
			std::vector<ALCTCandidate*> out_vec;

			/**********************
	 		* ALCT Candidate Setup
	 		**********************/

			for (int i=0; i<config.get_fifo_tbins()-config.get_drift_delay(); i++)
			{
				std::vector <ALCTCandidate*> temp_vec; 
				for (int j = 0; j<cvec.at(0)->get_maxWi(); j++)
				{
					ALCTCandidate * cand = new ALCTCandidate(j,1);
					temp_vec.push_back(cand);
				}
				end_vec.push_back(temp_vec); 
			}

			/**********************
	 		* RUNNING ALGORITHM
	 		**********************/

			trig_and_find(cvec, config, end_vec);
			ghostBuster(end_vec,config);
			extract_sort_cut(end_vec,out_vec);

			/**********************
	 		* MEMORY LEAK CLEANUP
	 		**********************/

			wipe(out_vec);
			wipe(cvec);

		}
		//alct_t_emu->Fill();
	}
	

	//outF->cd();
	//alct_t_emu->Write();

	//printf("Wrote to file: %s\n",outputfile.c_str());

	//auto t2 = std::chrono::high_resolution_clock::now();
	//cout << "Time elapsed: " << chrono::duration_cast<chrono::seconds>(t2-t1).count() << " s" << endl;
	return 0;
}