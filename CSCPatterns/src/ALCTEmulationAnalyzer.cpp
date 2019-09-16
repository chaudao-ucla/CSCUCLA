/*
 * ALCTEmulationAnalyzer.cpp
 *
 *  Created on: 12 September 2019
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

#include "../include/ALCTEmulationAnalyzer.h"

using namespace std; 

int main(int argc, char* argv[])
{
    ALCTEmulationAnalyzer p;
    return p.main(argc, argv);
}

int ALCTEmulationAnalyzer::run(string inputfile, int start, int end)
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

    /*TFile* outF = new TFile(outputfile.c_str(),"RECREATE");
	if(!outF)
	{
		printf("Failed to open output file: %s\n", outputfile.c_str());
		return -1;
	}

    TTree* alct_t_emu = new TTree("ALCTEmulationResults", "ALCTEmulationResults");*/

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

	int num_alct = 0; 
	int num_diff = 0; 
	int num_agree = 0; 

	int num_tot_alct=0; 
	int num_seg_unmatched = 0; 
	int num_seg_matched = 0;
	int num_alct_unmatched = 0; 
	int num_tot_seg=0;
	double mean_diff = 0.0;
	double std_diff = 0.0;
	int alct_match[4] = {0,0,0,0};
	int alct_unmatch[4] = {0,0,0,0};
	int track_match[3] = {0,0,0};
	int track_unmatch[3] = {0,0,0};
	bool flag_vec[alcts.size()];
	for (int i = 0; i< alcts.size(); i++)
	{
		flag_vec[i] = false; 
	}
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
			//extract(end_vec,out_vec);
			extract_sort_cut(end_vec,out_vec);

			num_tot_alct+=out_vec.size();

			/**********************
	 		* SEGMENT COMPARISON
	 		**********************/

		 	/*
			int cand_size = out_vec.size();

			for (int iter=0; iter<segments.size(); iter++)
			{
				int min_dist = 1000;
				int marker = - 1; 

				int chSid1 = CSCHelper::serialize(ST, RI, CH, EC);
				int chSid2 = chSid1;
				
				bool me11a	= ST == 1 && RI == 4;
				bool me11b	= ST == 1 && RI == 1;

				if (me11a || me11b)
				{
					if (me11a) chSid2 = CSCHelper::serialize(ST, 1, CH, EC);
					if (me11b) chSid2 = CSCHelper::serialize(ST, 4, CH, EC);
				}
				if (segments.mu_id->at(iter)==-1) continue; 
				if (chSid1!=segments.ch_id->at(iter) && chSid2!= segments.ch_id->at(iter)) continue;
				int pos_seg = segments.pos_y->at(iter)-1;
				if (pos_seg<0) continue; 
				num_tot_seg++;
				for (int j = 0; j<out_vec.size(); j++)
				{
					int temp_val = abs((int)(out_vec.at(j)->get_kwg()) - (int)pos_seg); 
					if (temp_val<min_dist)
					{
						min_dist = temp_val;
						marker = j;
					}
				}
				if (marker>-1)
				{
					double val = (double)((int)(out_vec.at(marker)->get_kwg()) - (int)pos_seg);
					if (min_dist>=2.0) cout << "Event = " << i << ", ST = " << ST << ", RI = " << RI << ", CH = " << CH << ", EC = " << EC << ", size = " << cand_size << ", min_dist = " << min_dist << endl << endl;
					mean_diff+= val;
					std_diff+= val*val;
					int temp = (int) (out_vec.at(marker)->get_quality());
					alct_match[temp]++;
					int temp2 = (int) (out_vec.at(marker)->get_tracknumber());
					track_match[temp2]++;
					out_vec.at(marker)->nix();
					out_vec.erase(out_vec.begin()+marker);
					num_seg_matched++;
					continue; 
				}
				num_seg_unmatched++;
			}

			for(int iter = 0; iter<out_vec.size(); iter++)
			{
				int temp = (int) (out_vec.at(iter)->get_quality());
				alct_unmatch[temp]++;
				int temp2 = (int) (out_vec.at(iter)->get_tracknumber());
				track_unmatch[temp2]++;
			}

			num_alct_unmatched+= out_vec.size();

			*/

			/**********************
	 		* EMULATION TESTER
	 		**********************/
			
			for (int i=0; i<alcts.size(); i++)
			{
				int chSid1 = CSCHelper::serialize(ST, RI, CH, EC);
				int chSid2 = chSid1;
				
				bool me11a	= ST == 1 && RI == 4;
				bool me11b	= ST == 1 && RI == 1;

				if (me11a || me11b)
				{
					if (me11a) chSid2 = CSCHelper::serialize(ST, 1, CH, EC);
					if (me11b) chSid2 = CSCHelper::serialize(ST, 4, CH, EC);
				}
				if (chSid1!=alcts.ch_id->at(i) && chSid2!= alcts.ch_id->at(i)) continue;
				for (int j = 0; j<out_vec.size(); j++)
				{
					ALCTCandidate* compare = out_vec[j];
					if (compare->get_first_bx() == alcts.BX->at(i)+5 && 
						compare->get_quality() >= alcts.quality->at(i) &&
						compare->get_kwg() == alcts.keyWG->at(i)
						)
					{
						num_agree++; 
						out_vec.at(j)->nix();
						out_vec.erase(out_vec.begin()+j);
					}
					else num_diff++; 
				}
				num_alct++;
			}
			
			/**********************
	 		* LCT TESTER (BEGIN)
	 		**********************/

			/*for (int i=0; i<lcts.size(); i++)
			{
				int chSid1 = CSCHelper::serialize(ST, RI, CH, EC);
				int chSid2 = chSid1;
				
				bool me11a	= ST == 1 && RI == 4;
				bool me11b	= ST == 1 && RI == 1;

				if (me11a || me11b)
				{
					if (me11a) chSid2 = CSCHelper::serialize(ST, 1, CH, EC);
					if (me11b) chSid2 = CSCHelper::serialize(ST, 4, CH, EC);
				}
				if (chSid1!=lcts.ch_id->at(i) && chSid2!= lcts.ch_id->at(i)) continue;
				num_lct++;
			}


			//if (num_alct!=out_vec.size()) cout << "Event = " << i << ", ST = " << ST << ", RI = " << RI << ", CH = " << CH << ", EC = " << EC << ", num_alct = " << num_alct << ", size = " << out_vec.size() << endl << endl;
			//if (num_alct==0 && candvec.size()) cout << "Event = " << i << ", ST = " << ST << ", RI = " << RI << ", CH = " << CH << ", EC = " << EC << ", num_alct = " << num_alct << ", size = " << candvec.size() << endl << endl;
			if (num_alct > 0 && num_lct == 0)
			{
				//if (out_vec.at(0)->get_first_bx()!=out_vec.at(1)->get_first_bx() && num_lct>=2)
				cout << "Event = " << i << ", ST = " << ST << ", RI = " << RI << ", CH = " << CH << ", EC = " << EC << ", num_alct = " << num_alct << ", size = " << out_vec.size() << endl << endl;
			}*/
			
			/*
			if (out_vec.size())
			{
				num_off+=out_vec.size(); 
				cout << "Event = " << i << ", ST = " << ST << ", RI = " << RI << ", CH = " << CH << ", EC = " << EC << endl << endl;
			}
			*/

			/**********************
	 		* LCT TESTER (END)
	 		**********************/

			/**********************
	 		* MEMORY LEAK (BEGIN)
	 		**********************/
			wipe(out_vec);
			wipe(cvec);
			/**********************
	 		* MEMORY LEAK (END)
	 		**********************/
		}
		//alct_t_emu->Fill();
	}
	
	mean_diff = mean_diff / (double) num_seg_matched;
	std_diff = sqrt((std_diff/(double)num_seg_matched) - mean_diff*mean_diff);

	cout<< "num_tot_alct = " <<  num_tot_alct << endl
		<< "num_tot_seg = " << num_tot_seg << endl
		<< "num_seg_matched = " << num_seg_matched << endl
		<< "num_seg_unmatched = " << num_seg_unmatched << endl
		<< "num_alct_unmatched = " << num_alct_unmatched << endl
		<< "mean_diff = " << mean_diff << endl
		<< "std_diff = " << std_diff << endl
		<< "matched q3 = " << alct_match[3] << endl
		<< "matched q2 = " << alct_match[2] << endl
		<< "matched q1 = " << alct_match[1] << endl
		<< "unmatched q3 = " << alct_unmatch[3] << endl
		<< "unmatched q2 = " << alct_unmatch[2] << endl
		<< "unmatched q1 = " << alct_unmatch[1] << endl
		<< "matched t1 = " << track_match[1] << endl
		<< "matched t2 = " << track_match[2] << endl
		<< "unmatched t1 = " << track_unmatch[1] << endl
		<< "unmatched t2 = " << track_unmatch[2] << endl
		<< endl << endl;
	

	cout << "num_alct = " << num_alct << endl
		 << "num_agree = " << num_agree << endl 
		 << "num_diff = " << num_diff << endl 
		 << endl << endl; 

	//outF->cd();
	//alct_t_emu->Write();

	//printf("Wrote to file: %s\n",outputfile.c_str());

	//auto t2 = std::chrono::high_resolution_clock::now();
	//cout << "Time elapsed: " << chrono::duration_cast<chrono::seconds>(t2-t1).count() << " s" << endl;
	return 0;
}