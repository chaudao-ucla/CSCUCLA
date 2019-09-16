#ifndef CSCUCLA_CSCDIGITUPLES_FILLCSCINFO_H
#define CSCUCLA_CSCDIGITUPLES_FILLCSCINFO_H


/* Original Author: Nick McColl
 *
 *   	    Author: William Nash
 *  Date Inherited: August 30th 2018
 *
 */

#include "FWCore/Framework/interface/Event.h"

#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>
#include <DataFormats/CSCDigi/interface/CSCStripDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCWireDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCComparatorDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCALCTDigiCollection.h>
#include "CSCUCLA/CSCDigiTuples/include/CSCHelper.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/EcalEndcapGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

//muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"

//segments
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "CSCInfo.h"


#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"

/* TODO:
 * - make TreeWriter, TreeReader class (inherit from TreeContainer?)
 * - make each info class (change name?) take as a parameter TreeWriter / TreeReader
 * 		-set branch address, or makes new branch on instantiation
 * 		-
 *
 * - fill should be overloaded:
 * 		- takes CMSSW variables with treewriter
 * 		- takes no argument with treereader, eats branch
 *
 *
 */

//for histogram bins, starting at 1
enum EVENT_CUTS {
	allEvents,
	hasSegments,
	hasCSCDigis,
	hasVertex,
	hasResonance,
	hasMuonInCSC,
	EVENT_SIZE
};

//for histogram bins, starting at 1
enum MUON_CUTS {
	allMuons,
	eventHasSegments,
	eventHasCSCDigis,
	eventHasVertex,
	isStandAlone,
	isGlobal,
	isInMassWindow,
	isOS,
	isOverPtThreshold,
	muonHasSegments,
	MUON_SIZE
};


typedef   unsigned char        size8 ; // 8 bit 0->255
typedef   unsigned short int   size16; //16 bit 0->65536
typedef   unsigned int         size  ; //32 bit 0->4294967296

class FillInfo;

/*
 * TODO: put all the histograms into this class, and write this on a tree write,
 * allows you to use multiple selections
 */

class SelectionHistograms {
public:
	SelectionHistograms(TreeContainer& t, const string& selection);
	~SelectionHistograms();
	void write();

	TH1F* h_eventCuts = 0;
	TH1F* h_muonCuts = 0;

	TH1F* h_osInvMass = 0;
	TH1F* h_ssInvMass = 0;
	TH1F* h_premassCutInvMass = 0;
	TH1F* h_nAllMuons = 0;
	TH1F* h_allMuonsPt = 0;
	TH1F* h_allMuonsEta = 0;
	TH1F* h_allMuonsPhi = 0;
	TH1F* h_allInvMass = 0;
	TH1F* h_selectedMuonsPt = 0;
	TH1F* h_selectedMuonsEta = 0;
	TH1F* h_selectedMuonsPhi = 0;
	TH1F* h_nSelectedMuons = 0;
	TH1F* h_nAllSegments = 0;
};

class TreeContainer {
public:
  TreeContainer(TString fileName, TString treeName, TString treeTitle, const string& selection):
  h(*this, selection){
	file = new TFile(fileName, "RECREATE");
	tree = new TTree(treeName,treeTitle);


  }
  void write() {
    file->cd();
    h.write();
    tree->Write();
    file->Close();
    delete file;
  }

  void fill();

  void reset();

  void addInfo(FillInfo* info){
	  if(info) infos.push_back(info);
  }

  TFile * file;
  TTree * tree;
  SelectionHistograms h;

private:
  std::vector<FillInfo*> infos;

};




class FillInfo {
public:
  FillInfo(const std::string& prfx, TreeContainer& tree) :
	  prefix(prfx),
	  fTree(&tree)
{
	  reset();
	  fTree->addInfo(this);
};
  virtual ~FillInfo() {};
  virtual void reset() {};
protected:
  //const TString prefix; //prefix before each tree branch associated with this object
  const string prefix;

  //Book single variable
  template<class T>
  void    book(const char *name, T& var, const char *type) { fTree->tree->Branch(
		  //TString(prefix).Append('_').Append(name).Data(),
		  (prefix+string("_")+string(name)).c_str(),
		 // TString(prefix).Append('_').Append(name).Data(),
		  &var,
		  (string(name)+string("/")+string(type)).c_str());}
		  //TString(name).Append("/").Append(type).Data()); }

  //Book vector
  template<class T>
  void    book(const char *name, std::vector<T>& varv)   { fTree->tree->Branch(
		  //TString(prefix).Append('_').Append(name).Data(),
		  (prefix+string("_")+string(name)).c_str(),
		  &varv); }

  TreeContainer * fTree;

};


class FillEventInfo : public CSCInfo::Event, public FillInfo {
public:

  FillEventInfo(TreeContainer& tree) :
	  Event(), FillInfo(name, tree)
{
    book(GET_VARIABLE_NAME(EventNumber),EventNumber,"l");
    book(GET_VARIABLE_NAME(RunNumber)  ,RunNumber  ,"l");
    book(GET_VARIABLE_NAME(LumiSection),LumiSection,"I");
    book(GET_VARIABLE_NAME(BXCrossing) ,BXCrossing ,"I");
    book(GET_VARIABLE_NAME(NSegmentsInEvent), NSegmentsInEvent, "I");
  }
  virtual ~FillEventInfo() {};

  virtual void reset(){
    EventNumber  = 0;
    RunNumber    = 0;
    LumiSection  = -1;
    BXCrossing   = -1;
    NSegmentsInEvent = 0;
  }

  void fill(const edm::Event& iEvent, unsigned int nSegments);

};

class FillGenParticleInfo : public CSCInfo::GenParticles, public FillInfo {
public:
	FillGenParticleInfo(TreeContainer& tree) :
		GenParticles(),
		FillInfo(name, tree)
{
		//might want to make a constructor for the object here...
		pdg_id = new std::vector<int>();
		pt = new std::vector<float>();
		eta = new std::vector<float>();
		phi = new std::vector<float>();
		q = new std::vector<int>();
		book(GET_VARIABLE_NAME(pdg_id), *pdg_id);
		book(GET_VARIABLE_NAME(pt), *pt);
		book(GET_VARIABLE_NAME(eta), *eta);
		book(GET_VARIABLE_NAME(phi), *phi);
		book(GET_VARIABLE_NAME(q), *q);
}
	virtual ~FillGenParticleInfo() {
		delete pdg_id;
		delete pt;
		delete eta;
		delete phi;
		delete q;
	}
	virtual void reset() {
		pdg_id->clear();
		pt->clear();
		eta->clear();
		phi->clear();
		q->clear();
	}

public:
	void fill(const vector<reco::GenParticle>& gen);
};

class FillSimHitsInfo : public CSCInfo::SimHits, public FillInfo {
public:
	FillSimHitsInfo(TreeContainer& tree) :
		SimHits(),
		FillInfo(name, tree)
{
		ch_id = new std::vector<int>();
		pdg_id = new std::vector<int>();
		layer = new std::vector<int>();
		energyLoss = new std::vector<float>();
		thetaAtEntry = new std::vector<float>();
		phiAtEntry = new std::vector<float>();
		pAtEntry = new std::vector<float>();
		book(GET_VARIABLE_NAME(ch_id), *ch_id);
		book(GET_VARIABLE_NAME(pdg_id), *pdg_id);
		book(GET_VARIABLE_NAME(layer), *layer);
		book(GET_VARIABLE_NAME(energyLoss), *energyLoss);
		book(GET_VARIABLE_NAME(thetaAtEntry), *thetaAtEntry);
		book(GET_VARIABLE_NAME(phiAtEntry), *phiAtEntry);
		book(GET_VARIABLE_NAME(pAtEntry), *pAtEntry);
}

	virtual ~FillSimHitsInfo() {
		delete ch_id;
		delete pdg_id;
		delete layer;
		delete energyLoss;
		delete thetaAtEntry;
		delete phiAtEntry;
		delete pAtEntry;
	}

	virtual void reset() {
		ch_id->clear();
		pdg_id->clear();
		layer->clear();
		energyLoss->clear();
		thetaAtEntry->clear();
		phiAtEntry->clear();
		pAtEntry->clear();
	}
public:
	void fill(const vector<PSimHit>& simhits);
};

class FillCaloHitsInfo : public CSCInfo::CaloHit, public FillInfo {
public:
	//pref = ecal, hcal, etc
	FillCaloHitsInfo(const string& pref, TreeContainer& tree) :
		CaloHit(pref),
		FillInfo(name,tree)
	{
		energyEM = new std::vector<float>();
		energyHad = new std::vector<float>();
		eta =  new std::vector<float>();
		phi = new std::vector<float>();
		book(GET_VARIABLE_NAME(energyEM), *energyEM);
		book(GET_VARIABLE_NAME(energyHad), *energyHad);
		book(GET_VARIABLE_NAME(eta), *eta);
		book(GET_VARIABLE_NAME(phi), *phi);
	}

	virtual ~FillCaloHitsInfo() {
		delete energyEM;
		delete energyHad;
		delete eta;
		delete phi;
	}

	virtual void reset() {
		energyEM->clear();
		energyHad->clear();
		eta->clear();
		phi->clear();
	}
public:
	virtual void fill(const vector<PCaloHit>& calHits, const EcalBarrelGeometry* theEcal) ;
	virtual void fill(const vector<PCaloHit>& calHits, const EcalEndcapGeometry* theEcal) ;
	virtual void fill(const vector<PCaloHit>& calHits, const EcalPreshowerGeometry* theEcal) ;
	virtual void fill(const vector<PCaloHit>& calHits, const HcalGeometry* theHcal) ;
};


class FillPFInfo : public CSCInfo::PFCandidate, public FillInfo {
public:
	FillPFInfo(TreeContainer& tree):
		PFCandidate(),
		FillInfo(name,tree)
	{
		pdg_id = new std::vector<int>();
		particleId = new std::vector<int>();
		eta = new std::vector<float>();
		phi = new std::vector<float>();
		ecalEnergy = new std::vector<float>();
		hcalEnergy = new std::vector<float>();
		h0Energy = new std::vector<float>();
		book(GET_VARIABLE_NAME(pdg_id), *pdg_id);
		book(GET_VARIABLE_NAME(particleId), *particleId);
		book(GET_VARIABLE_NAME(eta), *eta);
		book(GET_VARIABLE_NAME(phi), *phi);
		book(GET_VARIABLE_NAME(ecalEnergy), *ecalEnergy);
		book(GET_VARIABLE_NAME(hcalEnergy), *hcalEnergy);
		book(GET_VARIABLE_NAME(h0Energy), *h0Energy);
	}

	virtual ~FillPFInfo() {
		delete pdg_id;
		delete particleId;
		delete eta;
		delete phi;
		delete ecalEnergy;
		delete hcalEnergy;
		delete h0Energy;
	}

	virtual void reset() {
		pdg_id->clear();
		particleId->clear();
		eta->clear();
		phi->clear();
		ecalEnergy->clear();
		hcalEnergy->clear();
		h0Energy->clear();
	}
public:
	void fill(const vector<reco::PFCandidate>& pfCand);
};

class FillMuonInfo : public CSCInfo::Muons, public FillInfo {
public:
	FillMuonInfo(TreeContainer& tree) :
		Muons(),
		FillInfo(name,tree)
{

		//might want to make a constructor for the object here...
		pt = new std::vector<float>();
		eta = new std::vector<float>();
		phi = new std::vector<float>();
		q = new std::vector<int>();
		isGlobal = new std::vector<bool>();
		isTracker = new std::vector<bool>();
		book(GET_VARIABLE_NAME(pt), *pt);
		book(GET_VARIABLE_NAME(eta), *eta);
		book(GET_VARIABLE_NAME(phi), *phi);
		book(GET_VARIABLE_NAME(q), *q);
		book(GET_VARIABLE_NAME(isGlobal), *isGlobal);
		book(GET_VARIABLE_NAME(isTracker), *isTracker);
}
	virtual ~FillMuonInfo() {
		delete pt;
		delete eta;
		delete phi;
		delete q;
		delete isGlobal;
		delete isTracker;
	};

	virtual void reset() {
		pt->clear();
		eta->clear();
		phi->clear();
		q->clear();
		isGlobal->clear();
		isTracker->clear();
	}

public:

	void fill(const reco::MuonCollection& muons);
	void fill(const reco::Muon& muon);
};


class FillSegmentInfo : public CSCInfo::Segments, public FillInfo {
public:

  FillSegmentInfo(TreeContainer& tree) :
	  Segments(),
	  FillInfo(name,tree) {
	mu_id = new vector<int>();
	ch_id = new vector<int>();
	pos_x = new vector<float>();
	pos_y = new vector<int>();
	dxdz = new vector <float>();
	dydz = new vector<float>();
	chisq = new vector<float>();
	nHits = new vector<size8>();
    book(GET_VARIABLE_NAME(mu_id)     ,*mu_id      );
    book(GET_VARIABLE_NAME(ch_id)     ,*ch_id      );
    book(GET_VARIABLE_NAME(pos_x)	  ,*pos_x      );
    book(GET_VARIABLE_NAME(pos_y)     ,*pos_y      );
    book(GET_VARIABLE_NAME(dxdz)      ,*dxdz       );
    book(GET_VARIABLE_NAME(dydz)      ,*dydz       );
    book(GET_VARIABLE_NAME(chisq)     ,*chisq      );
    book(GET_VARIABLE_NAME(nHits)     ,*nHits      );

  }
  virtual ~FillSegmentInfo() {
	  delete mu_id;
	  delete ch_id;
	  delete pos_x;
	  delete pos_y;
	  delete dxdz;
	  delete dydz;
	  delete chisq;
	  delete nHits;
  };

  virtual void reset(){
    mu_id       ->clear();
    ch_id       ->clear();
    pos_x       ->clear();
    pos_y       ->clear();
    dxdz        ->clear();
    dydz        ->clear();
    chisq       ->clear();
    nHits       ->clear();
  }

  //void fill(const CSCSegmentCollection& segments, const CSCRecHit2DCollection * recHits = 0);
  void fill(std::vector<const CSCSegment*>& segments, const CSCGeometry* theCSC, int mu_index = -1);
  void fill(const CSCSegment& segment, const CSCGeometry* theCSC, int mu_index);
  size16 findRecHitIdx(const CSCRecHit2D& hit, const CSCRecHit2DCollection* allRecHits);

};

class FillRecHitInfo : public CSCInfo::RecHits, public FillInfo {
public:

	FillRecHitInfo(TreeContainer& tree) :
			RecHits(),
			FillInfo(name, tree) {
		mu_id = new vector<int>();
		ch_id = new vector<int>();
		lay = new vector<size8>();
		pos_x = new vector<float>();
		pos_y = new vector<float>();
		e = new vector<float>();
		max_adc = new vector<float>();
		book(GET_VARIABLE_NAME(mu_id), *mu_id);
		book(GET_VARIABLE_NAME(ch_id), *ch_id);
		book(GET_VARIABLE_NAME(lay), *lay);
		book(GET_VARIABLE_NAME(pos_x), *pos_x);
		book(GET_VARIABLE_NAME(pos_y), *pos_y);
		book(GET_VARIABLE_NAME(e), *e);
		book(GET_VARIABLE_NAME(max_adc), *max_adc);

	}
	virtual ~FillRecHitInfo() {
		delete mu_id;
		delete ch_id;
		delete lay;
		delete pos_x;
		delete pos_y;
		delete e;
		delete max_adc;
	}
	;

	virtual void reset() {
		mu_id->clear();
		ch_id->clear();
		lay->clear();
		pos_x->clear();
		pos_y->clear();
		e->clear();
		max_adc->clear();
	}

  void fill(const std::vector<CSCRecHit2D>& recHits, int mu_index);

};

class FillLCTInfo: public CSCInfo::LCTs, public FillInfo {
public:

	FillLCTInfo(TreeContainer& tree) :
			LCTs(),
			FillInfo(name, tree) {
		ch_id = new vector<size16>();
		quality = new vector<size8>();
		pattern = new vector<size8>();
		bend = new vector<size8>();
		keyWireGroup = new vector<size8>();
		keyHalfStrip = new vector<size8>();
		bunchCross = new vector<size8>();
		book(GET_VARIABLE_NAME(ch_id), *ch_id);
		book(GET_VARIABLE_NAME(quality), *quality);
		book(GET_VARIABLE_NAME(pattern), *pattern);
		book(GET_VARIABLE_NAME(bend), *bend);
		book(GET_VARIABLE_NAME(keyWireGroup), *keyWireGroup);
		book(GET_VARIABLE_NAME(keyHalfStrip), *keyHalfStrip);
		book(GET_VARIABLE_NAME(bunchCross), *bunchCross);

	}
	virtual ~FillLCTInfo() {
		delete ch_id;
		delete quality;
		delete pattern;
		delete bend;
		delete keyWireGroup;
		delete keyHalfStrip;
		delete bunchCross;
	};

	virtual void reset() {
		ch_id->clear();
		quality->clear();
		pattern->clear();
		bend->clear();
		keyWireGroup->clear();
		keyHalfStrip->clear();
		bunchCross->clear();
	}

	void fill(const CSCCorrelatedLCTDigiCollection& lcts);

};



class FillCLCTInfo: public CSCInfo::CLCTs, public FillInfo {
public:

	FillCLCTInfo(TreeContainer& tree) :
			CLCTs(),
			FillInfo(name, tree) {
		ch_id = new vector<size16>();
		isValid = new vector<size8>();
		quality = new vector<size16>();
		pattern = new vector<size8>();
		stripType = new vector<size8>();
		bend = new vector<size8>();
		halfStrip = new vector<size8>();
		CFEB = new vector<size8>();
		BX = new vector<size8>();
		trkNumber = new vector<size8>();
		keyStrip = new vector<size8>(); //this is in half strips for some reason, following cmssw convention
		book(GET_VARIABLE_NAME(ch_id), *ch_id);
		book(GET_VARIABLE_NAME(isValid), *isValid);
		book(GET_VARIABLE_NAME(quality), *quality);
		book(GET_VARIABLE_NAME(pattern), *pattern);
		book(GET_VARIABLE_NAME(stripType), *stripType);
		book(GET_VARIABLE_NAME(bend), *bend);
		book(GET_VARIABLE_NAME(halfStrip), *halfStrip);
		book(GET_VARIABLE_NAME(CFEB), *CFEB);
		book(GET_VARIABLE_NAME(BX), *BX);
		book(GET_VARIABLE_NAME(trkNumber), *trkNumber);
		book(GET_VARIABLE_NAME(keyStrip), *keyStrip);

	}
	virtual ~FillCLCTInfo() {
		delete ch_id;
		delete isValid;
		delete quality;
		delete pattern;
		delete stripType;
		delete bend;
		delete halfStrip;
		delete CFEB;
		delete BX;
		delete trkNumber;
		delete keyStrip;
	}

	virtual void reset() {
		ch_id->clear();
		isValid->clear();
		quality->clear();
		pattern->clear();
		stripType->clear();
		bend->clear();
		halfStrip->clear();
		CFEB->clear();
		BX->clear();
		trkNumber->clear();
		keyStrip->clear();
	}

	void fill(const CSCCLCTDigiCollection& clcts);

};

class FillALCTInfo : public CSCInfo::ALCTs, public FillInfo 
{
	public: 
		FillALCTInfo(TreeContainer& tree) :
			ALCTs(),
			FillInfo(name,tree)
		{
			ch_id = new std::vector<size16>();
			isValid = new std::vector<size8>();
			quality = new std::vector<size8>();
			accelerator = new std::vector<size8>();
			collisionB = new std::vector<size8>();
			keyWG = new std::vector<size8>();
			BX = new std::vector<size8>();
			trkNumber = new std::vector<size8>();
			fullBX = new std::vector<size16>();

			book(GET_VARIABLE_NAME(ch_id),*ch_id);
			book(GET_VARIABLE_NAME(isValid),*isValid);
			book(GET_VARIABLE_NAME(quality),*quality);
			book(GET_VARIABLE_NAME(accelerator),*accelerator);
			book(GET_VARIABLE_NAME(collisionB),*collisionB);
			book(GET_VARIABLE_NAME(keyWG),*keyWG);
			book(GET_VARIABLE_NAME(BX),*BX);
			book(GET_VARIABLE_NAME(trkNumber),*trkNumber);
			book(GET_VARIABLE_NAME(fullBX),*fullBX);
		}

		virtual ~FillALCTInfo()
		{
			delete ch_id;
			delete isValid;
			delete quality;
			delete accelerator; 
			delete collisionB;
			delete keyWG;
			delete BX;
			delete trkNumber;
			delete fullBX;
		}

		virtual void reset()
		{
			ch_id->clear();
			isValid->clear();
			quality->clear();
			accelerator->clear(); 
			collisionB->clear();
			keyWG->clear();
			BX->clear();
			trkNumber->clear();
			fullBX->clear();
		}

		void fill(const CSCALCTDigiCollection &alcts);
};

class FillWireInfo : public CSCInfo::Wires, public FillInfo 
{
	public: 
		FillWireInfo(TreeContainer& tree) :
			Wires(),
			FillInfo(name,tree)
		{
			ch_id = new std::vector<size16>();
			group = new std::vector<int>();
			lay = new std::vector<size8>();
			timeBin = new std::vector<int>();
			BX = new std::vector<int>();
			timeBinWord = new std::vector<uint32_t>();
			//time_bins_on = new std::vector<std::vector<int>>();

			book(GET_VARIABLE_NAME(ch_id),*ch_id);
			book(GET_VARIABLE_NAME(group),*group);
			book(GET_VARIABLE_NAME(lay),*lay);
			book(GET_VARIABLE_NAME(timeBin),*timeBin);
			book(GET_VARIABLE_NAME(timeBinWord),*timeBinWord);
			book(GET_VARIABLE_NAME(BX),*BX);
			//book(GET_VARIABLE_NAME(time_bins_on),*time_bins_on);
		}

		virtual ~FillWireInfo()
		{
			delete ch_id;
			delete group;
			delete lay;
			delete timeBin;
			delete timeBinWord;
			delete BX;
			//delete time_bins_on;
		}

		virtual void reset()
		{
			ch_id->clear();
			group->clear();
			lay->clear();
			timeBin->clear();
			timeBinWord->clear();
			BX->clear();
			//time_bins_on->clear();
		}

		void fill(const CSCWireDigiCollection &wires);
};

class FillStripInfo : public CSCInfo::Strips, public FillInfo
{
	public:
		FillStripInfo(TreeContainer &tree):
			Strips(),
			FillInfo(name,tree)
		{
			ch_id = new std::vector<size16>();
			lay = new std::vector<size8>();
			num = new std::vector<int>();
//			ADC = new std::vector<std::vector<int>>();
//			L1APhase = new std::vector<std::vector<int>>();
//			ADCOverflow = new std::vector<std::vector<uint16_t>>();
//			OverlappedSample = new std::vector<std::vector<uint16_t>>();
//			Errorstat = new std::vector<std::vector<uint16_t>>();

			book(GET_VARIABLE_NAME(ch_id), *ch_id);
			book(GET_VARIABLE_NAME(lay), *lay);
			book(GET_VARIABLE_NAME(num), *num);
//			book(GET_VARIABLE_NAME(ADC), *ADC);
//			book(GET_VARIABLE_NAME(L1APhase), *L1APhase);
//			book(GET_VARIABLE_NAME(ADCOverflow), *ADCOverflow);
//			book(GET_VARIABLE_NAME(OverlappedSample), *OverlappedSample);
//			book(GET_VARIABLE_NAME(Errorstat), *Errorstat);
		}

		virtual ~FillStripInfo()
		{
			delete ch_id;
			delete lay;
			delete num;
//			delete ADC;
//			delete L1APhase;
//			delete ADCOverflow;
//			delete OverlappedSample;
//			delete Errorstat;
		}

		virtual void reset()
		{
			ch_id->clear();
			lay->clear();
			num->clear();
//			ADC->clear();
//			L1APhase->clear();
//			ADCOverflow->clear();
//			OverlappedSample->clear();
//			Errorstat->clear();
		}

		void fill(const CSCStripDigiCollection &strips);
};

class FillCompInfo : public CSCInfo::Comparators, public FillInfo {
public:

  FillCompInfo(TreeContainer& tree) :
	  Comparators(),
	  FillInfo(name, tree) {
	  ch_id = new vector<int>();
	  lay = new vector<size8>();
	  strip = new vector<size8>();
	  halfStrip = new vector<size8>();
	  bestTime = new vector<size8>();
	  nTimeOn = new vector<size8>();

		book(GET_VARIABLE_NAME(ch_id), *ch_id);
		book(GET_VARIABLE_NAME(lay), *lay);
		book(GET_VARIABLE_NAME(strip), *strip);
		book(GET_VARIABLE_NAME(halfStrip), *halfStrip);
		book(GET_VARIABLE_NAME(bestTime), *bestTime);
		book(GET_VARIABLE_NAME(nTimeOn), *nTimeOn);
  }

  virtual ~FillCompInfo() {
	  delete ch_id;
	  delete lay;
	  delete strip;
	  delete halfStrip;
	  delete bestTime;
	  delete nTimeOn;
  };

  virtual void reset(){
	  ch_id->clear();
	  lay->clear();
	  strip->clear();
	  halfStrip->clear();
	  bestTime->clear();
	  nTimeOn->clear();
  }

  void fill(const CSCComparatorDigiCollection& strips);

};




/*
class FillStripInfo : public FillInfo {
public:

  FillStripInfo(TreeContainer& tree) :FillInfo(tree) {
    book("strip_id"           ,strip_id           );
    book("strip_lay"          ,strip_lay          );
    book("strip_number"       ,strip_number      );



  }
  virtual ~FillStripInfo() {};

private:
  std::vector<size16>   strip_id     ;
  std::vector<size8>    strip_lay    ;
  std::vector<size8>    strip_number ;


  virtual void reset(){
    strip_id     .clear();
    strip_lay    .clear();
    strip_number .clear();
  }

  public:

  void fill(const CSCStripDigiCollection& strips);

};

class FillCompInfo : public FillInfo {
public:

  FillCompInfo(TreeContainer& tree) :FillInfo(tree) {
    book("comp_id"      ,comp_id          );
    book("comp_lay"     ,comp_lay         );
    book("comp_strip"   ,comp_strip       );
    book("comp_comp"    ,comp_comp        );



  }
  virtual ~FillCompInfo() {};

private:
  std::vector<size16> comp_id    ;
  std::vector<size8>  comp_lay   ;
  std::vector<size8>  comp_strip ;
  std::vector<size8>  comp_comp  ;


  virtual void reset(){
    comp_id    .clear();
    comp_lay   .clear();
    comp_strip .clear();
    comp_comp  .clear();
  }

  public:

  void fill(const CSCComparatorDigiCollection& strips);

};


class FillWireInfo : public FillInfo {
public:

  FillWireInfo(TreeContainer& tree) :FillInfo(tree) {
    book("wire_id"     ,wire_id         );
    book("wire_lay"    ,wire_lay        );
    book("wire_grp"    ,wire_grp        );
    book("wire_time"   ,wire_time       );
    book("wire_bx"     ,wire_bx         );



  }
  virtual ~FillWireInfo() {};

private:
  std::vector<size16>     wire_id  ;
  std::vector<size8>      wire_lay ;
  std::vector<size8>      wire_grp ;
  std::vector<size8>      wire_time;
  std::vector<int>        wire_bx;


  virtual void reset(){
    wire_id   .clear();
    wire_lay  .clear();
    wire_grp  .clear();
    wire_time .clear();
    wire_bx   .clear();
  }

  public:

  void fill(const CSCWireDigiCollection& wires);

};





class FillALCTInfo : public FillInfo {
public:

  FillALCTInfo(TreeContainer& tree) :FillInfo(tree) {

    book("alct_id"         , alct_id        );
    book("alct_isvalid"    , alct_isvalid   );
    book("alct_quality"    , alct_quality   );
    book("alct_accel"      , alct_accel     );
    book("alct_collB"      , alct_collB     );
    book("alct_wireGroup"  , alct_wireGroup );
    book("alct_BX"         , alct_BX        );
    book("alct_trkNumber"  , alct_trkNumber );



  }
  virtual ~FillALCTInfo() {};

private:
   std::vector<size16> alct_id          ;
   std::vector<size8>  alct_isvalid     ;
   std::vector<size8>  alct_quality     ;
   std::vector<size8>  alct_accel       ;
   std::vector<size8>  alct_collB       ;
   std::vector<size8>  alct_wireGroup   ;
   std::vector<size8>  alct_BX          ;
   std::vector<size8>  alct_trkNumber   ;


  virtual void reset(){
    alct_id         .clear();
    alct_isvalid    .clear();
    alct_quality    .clear();
    alct_accel      .clear();
    alct_collB      .clear();
    alct_wireGroup  .clear();
    alct_BX         .clear();
    alct_trkNumber  .clear();
  }

  public:

  void fill(const CSCALCTDigiCollection& alcts);

};

*/

#endif /*CSCUCLA_CSCDIGITUPLES_FILLCSCINFO_H*/
