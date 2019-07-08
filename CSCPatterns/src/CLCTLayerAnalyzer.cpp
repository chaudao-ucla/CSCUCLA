/*
 * CLCTLayerAnalyzer.cpp
 *
 *  Created on: Oct 18, 2018
 *      Author: wnash
 */


#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include <string>
#include <stdio.h>
#include <algorithm>
#include <time.h>

using namespace std;

#include "../include/CSCConstants.h"
#include "../include/CSCClasses.h"
#include "../include/CSCHelperFunctions.h"
#include "../include/LUTClasses.h"
#include "../include/CSCInfo.h"
#include "../include/CSCHelper.h"

#include "../include/CLCTLayerAnalyzer.h"

int main(int argc, char* argv[]){
	CLCTLayerAnalyzer p;
	return p.main(argc,argv);
}


//looks at a map of histograms, if it contains the key, fills the corresponding histogram with "histValue"
void CLCTLayerAnalyzer::fillHist(map<unsigned int, TH1F*> hists, unsigned int key, float histValue){
	//look to see if we care about this chamber
	auto it = hists.find(key);
	if(it != hists.end()){
		//fill the correct histogram
		it->second->Fill(histValue);
	}
}

map<unsigned int, TH1F*> CLCTLayerAnalyzer::makeHistPermutation(string name, string title, unsigned int bins, unsigned int low, unsigned int high){
	TH1F* h_me11a_plus = new TH1F((name+"_me_p11a_11").c_str(),("me_p11a_11 " +title).c_str(),bins,low,high);
	TH1F* h_me11b_plus = new TH1F((name+"_me_p11b_11").c_str(),("me_p11b_11" +title).c_str(),bins,low,high);
	TH1F* h_me11a_minus = new TH1F((name+"_me_m11a_11").c_str(),("me_m11a_11 "+title).c_str(),bins,low,high);
	TH1F* h_me11b_minus = new TH1F((name+"_me_m11b_11").c_str(),("me_m11b_11" +title).c_str(),bins,low,high);

	unsigned int me_p11a11 = CSCHelper::serialize(1,4,11,1);
	unsigned int me_p11b11 = CSCHelper::serialize(1,1,11,1);
	unsigned int me_m11a11 = CSCHelper::serialize(1,4,11,2);
	unsigned int me_m11b11 = CSCHelper::serialize(1,1,11,2);

	map<unsigned int, TH1F*> hists;
	hists[me_p11a11] = h_me11a_plus;
	hists[me_p11b11] = h_me11b_plus;
	hists[me_m11a11] = h_me11a_minus;
	hists[me_m11b11] = h_me11b_minus;
	return hists;
}


int CLCTLayerAnalyzer::run(string inputfile, string outputfile, int start, int end) {


	cout << "Running over file: " << inputfile << endl;


	TFile* f = TFile::Open(inputfile.c_str());
	if(!f) throw "Can't open file";

	TTree* t =  (TTree*)f->Get("CSCDigiTree");
	if(!t) throw "Can't find tree";


	//
	// SET INPUT BRANCHES
	//

	CSCInfo::Event evt(t);
	CSCInfo::Muons muons(t);
	CSCInfo::Segments segments(t);
	CSCInfo::RecHits recHits(t);
	CSCInfo::LCTs lcts(t);
	CSCInfo::CLCTs clcts(t);
	CSCInfo::Comparators comparators(t);



	TFile * outF = new TFile(outputfile.c_str(),"RECREATE");
	if(!outF){
		printf("Failed to open output file: %s\n", outputfile.c_str());
		return -1;
	}


	map<unsigned int,TH1F*> allClctLayerCounts = makeHistPermutation("h_allClctLayerCount", "h_allClctLayerCount;Layer Count; CLCTs",7,0,7);
	map<unsigned int,TH1F*> clctLayerCounts = makeHistPermutation("h_clctLayerCount", "h_clctLayerCount;Layer Count; Matched CLCTs",7,0,7);
	map<unsigned int,TH1F*> unmatchedClctLayerCounts = makeHistPermutation("h_unmatchedClctLayerCount", "h_unmatchedClctLayerCount;Layer Count; Unmatched CLCTs",7,0,7);


	//all segments in the chamber with the clct layer threshold changed
	map<unsigned int,TH1F*> clctEff_den = makeHistPermutation("h_clctEff_den", "h_clctEff_den; Pt [GeV]; Efficiency", 16, 20,100);
	//if the chamber associated with the segment has any clct in it
	map<unsigned int,TH1F*> clctEff_hasClct = makeHistPermutation("h_clctEff_hasClct", "h_clctEff_hasClct; Pt [GeV]; Efficiency", 16, 20,100);
	//if the chamber has a 3 layer clct in it
	map<unsigned int,TH1F*>clctEff_has3LayClct = makeHistPermutation("h_clctEff_has3layClct", "h_clctEff_has3LayClct; Pt [GeV]; Efficiency", 16, 20,100);
	//if the 3 layer clct is the closest to the segment
	map<unsigned int,TH1F*> clctEff_3LayClct = makeHistPermutation("h_clctEff_3LayClct", "h_clctEff_3LayClct; Pt [GeV]; Efficiency", 16, 20,100);


	enum CLCT_EFF_CUTS {
		nSegments,
		hasClct,
		has3LayClct,
		match3LayClct
	};
	map<unsigned int,TH1F*> clctEff_cuts = makeHistPermutation("h_clctEff_cuts", "h_clctEff_3LayClct;; Segments", 4, 0,4);
	for(auto it : clctEff_cuts){
		it.second->GetXaxis()->SetBinLabel(nSegments+1, "nSegments");
		it.second->GetXaxis()->SetBinLabel(hasClct+1, "hasClct");
		it.second->GetXaxis()->SetBinLabel(has3LayClct+1, "has3LayClct");
		it.second->GetXaxis()->SetBinLabel(match3LayClct+1, "match3LayClct");
	}


	map<unsigned int, TH1F*> lctQuality = makeHistPermutation("h_lctQuality","h_lctQuality; Quality;LCTs",16,0,16);
	map<unsigned int, TH1F*> allLctQuality = makeHistPermutation("h_allLctQuality","h_allLctQuality; Quality;LCTs",16,0,16);

	map<unsigned int, TH1F*> matchPt = makeHistPermutation("h_matchedPt", "h_matchedPt;Pt [GeV]; CLCTs", 18,0,90);
	map<unsigned int, TH1F*> matchPt_3Lay = makeHistPermutation("h_matchedPt_3Lay", "h_matchedPt_3Lay;Pt [GeV]; CLCTs", 18,0,90);


	//int EC = 0; // 1-2
	int ST = 0; // 1-4
	int RI = 0; // 1-4
	//int CH = 0;
	float segmentX = 0;

	if(end > t->GetEntries() || end < 0) end = t->GetEntries();

	printf("Starting Event = %i, Ending Event = %i\n", start, end);
	//t->SetImplicitMT(true);


	unsigned int threeLayerCounterME11A = 0;

	for(int i = start; i < end; i++) {
		if(!(i%10000)) printf("%3.2f%% Done --- Processed %u Events\n", 100.*(i-start)/(end-start), i-start);

		t->GetEntry(i);
		//
		// Oct. 29 Skipping all events past new firmware update
		//
		//
		// Nov. 6 CHANGE ME BACK, NOW LOOKING AT LATER ERA!
		//
		/* First 3-layer firmware installation era on ME+1/1/11. Does not include min-CLCT-separation change (10 -> 5)
		 * installed on September 12
		 */
		if(evt. RunNumber < 321710 || evt.RunNumber > 323362) continue; //correct
		/* Era after min-separation change (10 -> 5), also includes 3 layer firmware change
		 */
		//if(evt.RunNumber <= 323362) continue;


		vector< unsigned int> matchedCLCTs;
		vector<unsigned int> matchedLCTs;

		//iterate through segments
		for(unsigned int thisSeg = 0; thisSeg < segments.size(); thisSeg++){
			int segId = segments.ch_id->at(thisSeg);
			CSCHelper::ChamberId c = CSCHelper::unserialize(segId);

			//EC = c.endcap;
			ST = c.station;
			RI = c.ring;
			//CH = c.chamber;


			segmentX = segments.pos_x->at(thisSeg); //strips


			// IGNORE SEGMENTS AT THE EDGES OF THE CHAMBERS
			if(CSCHelper::segmentIsOnEdgeOfChamber(segmentX, ST,RI)) continue;

			bool me11a = (ST == 1 && RI == 4);
			bool me11b = (ST == 1 && RI == 1);
			if(!(me11b || me11a)) continue;

			//Selecting only CLCTs in these chambers
			if(segments.mu_id->at(thisSeg) == -1) continue; //skip segments without a muon
			float Pt = muons.pt->at(segments.mu_id->at(thisSeg));

			/*
			 * Cutting on Pt!
			 */
			if (Pt < 25) continue;

			fillHist(clctEff_den, segId, Pt);
			fillHist(clctEff_cuts,segId,nSegments);

			bool foundCLCT = false; //if we found a clct in the same chamber as the segment
			bool found3LayCLCT = false;
			bool matched3LayCLCT = false;

			int closestCLCTtoSegmentIndex = -1;
			float minDistanceSegmentToClosestCLCT = 1e5;
			for(unsigned int iclct =0; iclct < clcts.size(); iclct++){
				int thisClctId = clcts.ch_id->at(iclct);
				if(thisClctId != segId) continue;
				foundCLCT = true;
				unsigned int qual = clcts.quality->at(iclct);
				if(qual == 3) found3LayCLCT = true;


				if(std::find(matchedCLCTs.begin(), matchedCLCTs.end(), iclct) != matchedCLCTs.end()) continue;
				float clctStripPos = clcts.halfStrip->at(iclct) / 2. + 16*clcts.CFEB->at(iclct);
				if(me11a) clctStripPos -= 16*4;
				if(abs(clctStripPos - segmentX) < minDistanceSegmentToClosestCLCT) {
					minDistanceSegmentToClosestCLCT = abs(clctStripPos - segmentX);
					closestCLCTtoSegmentIndex = iclct;
				}
			}
			if(closestCLCTtoSegmentIndex != -1){ //if we found one
				matchedCLCTs.push_back((unsigned int)closestCLCTtoSegmentIndex);
				unsigned int qual = clcts.quality->at(closestCLCTtoSegmentIndex);

				if(qual == 3) {
					cout << "threeLayerCount = " << ++threeLayerCounterME11A << endl;
					//printChamber(theseCompHits);
				}

				//look to see if we care about this chamber
				fillHist(clctLayerCounts,segId, qual);

				fillHist(matchPt,segId, Pt);
				if(qual == 3) fillHist(matchPt_3Lay,segId,Pt);

				//if(qual == 3) printChamber(theseCompHits);
				if(qual == 3) matched3LayCLCT = true;

				//closest index
			}else {
				fillHist(clctLayerCounts,segId, 0); //set a quality 0 if we find no associated clct

			}
			if(foundCLCT) {
				fillHist(clctEff_hasClct,segId,Pt);
				fillHist(clctEff_cuts,segId,hasClct);

			}
			if(found3LayCLCT) {
				fillHist(clctEff_has3LayClct,segId,Pt);
				fillHist(clctEff_cuts,segId, has3LayClct);

			}
			if(matched3LayCLCT) {
				fillHist(clctEff_3LayClct,segId,Pt);
				fillHist(clctEff_cuts,segId, match3LayClct);
			}


			//MATCH LCTS
			int closestLCTtoSegmentIndex = -1;
			float minDistanceSegmentToClosestLCT = 1e5;
			for(unsigned int ilct =0; ilct < lcts.size(); ilct++){
				int thisLctId = lcts.ch_id->at(ilct);
				if(thisLctId != segId) continue;
				if(std::find(matchedLCTs.begin(), matchedLCTs.end(),  ilct) != matchedLCTs.end()) continue;
				float lctStripPos = lcts.keyHalfStrip->at(ilct)/2.;
				if(me11a) lctStripPos -= 16*4;
				//if(sameChamberCount > 1)cout << "segPos = "<< segmentX <<" lctPos = " << lctStripPos << endl;
				if(abs(lctStripPos - segmentX) < minDistanceSegmentToClosestLCT) {
					minDistanceSegmentToClosestLCT = abs(lctStripPos - segmentX);
					closestLCTtoSegmentIndex = ilct;
				}
			}
			if(closestLCTtoSegmentIndex != -1){ //if we found one
				matchedLCTs.push_back((unsigned int)closestLCTtoSegmentIndex);
				unsigned int qual = lcts.quality->at(closestLCTtoSegmentIndex);

				//look to see if we care about this chamber
				fillHist(lctQuality,segId, qual);
			} else {
				fillHist(lctQuality,segId,0); //call unmatched "quality 0"
			}
		}//segments


		//now look through the clcts that were not matched to a segment
		for(unsigned int iclct =0; iclct < clcts.size(); iclct++){
			int thisClctId = clcts.ch_id->at(iclct);
			unsigned int qual = clcts.quality->at(iclct);

			fillHist(allClctLayerCounts,thisClctId,qual);

			if(std::find(matchedCLCTs.begin(), matchedCLCTs.end(), iclct) != matchedCLCTs.end()) continue;

			fillHist(unmatchedClctLayerCounts,thisClctId,qual);

		}

		for(unsigned int ilct=0; ilct < lcts.size(); ilct++){
			int thisLctId = lcts.ch_id->at(ilct);
			int qual = lcts.quality->at(ilct);
			fillHist(allLctQuality, thisLctId, qual);
		}
	}


	for(auto hist : allClctLayerCounts) hist.second->Write();
	for(auto hist : clctLayerCounts) hist.second->Write();
	for(auto hist : unmatchedClctLayerCounts) hist.second->Write();
	for(auto hist : clctEff_den) hist.second->Write();
	for(auto hist : clctEff_hasClct) hist.second->Write();
	for(auto hist : clctEff_has3LayClct) hist.second->Write();
	for(auto hist : clctEff_3LayClct) hist.second->Write();
	for(auto hist : clctEff_cuts) hist.second->Write();
	for(auto hist : lctQuality) hist.second->Write();
	for(auto hist : allLctQuality) hist.second->Write();
	for(auto hist : matchPt) hist.second->Write();
	for(auto hist : matchPt_3Lay) hist.second->Write();

	outF->Close();

	cout << "Wrote to file: " << outputfile << endl;

	return 0;
}


