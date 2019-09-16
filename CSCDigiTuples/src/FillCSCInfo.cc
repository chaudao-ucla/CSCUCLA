#include "CSCUCLA/CSCDigiTuples/include/FillCSCInfo.h"
#include "CSCUCLA/CSCDigiTuples/include/CSCHelper.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"


void TreeContainer::fill()
{
	  tree->Fill();
}

void TreeContainer::reset()
{
	for(auto& info: infos) info->reset();
}

SelectionHistograms::SelectionHistograms(TreeContainer& t, const string& selection){

	h_eventCuts = new TH1F("h_eventCuts", "h_eventCuts; ; Events", EVENT_CUTS::EVENT_SIZE, 0, EVENT_CUTS::EVENT_SIZE);
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::allEvents+1, "allEvents");
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasSegments+1, "eventHasSegments");
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasCSCDigis+1, "eventHasCSCDigis");
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasVertex+1, "eventHasVertex");
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasResonance+1, "eventHasResonance");
	h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasMuonInCSC+1, "eventHasMuonInCSC");


	h_muonCuts = new TH1F("h_muonCuts", "h_muonCuts; ; Muons", MUON_CUTS::MUON_SIZE, 0, MUON_CUTS::MUON_SIZE);
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::allMuons+1, "allMuons");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::eventHasSegments+1, "eventHasSegments");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::eventHasCSCDigis+1, "eventHasCSCDigis");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::eventHasVertex+1, "eventHasVertex");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::isStandAlone+1, "isStandalone");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::isGlobal+1, "isGlobal");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::isInMassWindow+1, "isInMassWindow");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::isOS+1, "isOS");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::isOverPtThreshold+1, "isOverPtThreshold");
	h_muonCuts->GetXaxis()->SetBinLabel(MUON_CUTS::muonHasSegments+1, "muonHasSegment");

	h_allMuonsEta = new TH1F("h_allMuonsEta", "h_allMuonsEta; #eta; Muons", 100, -3,3);
	h_allMuonsPhi = new TH1F("h_allMuonsPhi", "h_allMuonsPhi; #phi; Muons", 100, -3.2,3.2);
	h_allInvMass = new TH1F("h_allInvMass", "h_allInvMass; Mass [GeV], Dimuons", 130, 0, 130);

	h_selectedMuonsEta = new TH1F("h_selectedMuonsEta", "h_selectedMuonsEta; #eta; Muons", 100, -3,3);
	h_selectedMuonsPhi = new TH1F("h_selectedMuonsPhi", "h_selectedMuonsPhi; #phi; Muons", 100, -3.2,3.2);
	h_nAllMuons = new TH1F("h_nAllMuons", "h_nAllMuons; Muons; Events", 20,0,20);
	h_nSelectedMuons = new TH1F("h_nSelectedMuons", "h_nSelectedMuons; Muons; Events", 20, 0, 20);
	h_nAllSegments = new TH1F("h_nAllSegments", "h_nAllSegments; Segments; Count", 25, 0,25);


	int invMassBins = 100;
	if (selection == "ZeroBias"){
		cout <<  "--- Tree instantiated as minBias sample --- " << endl;
		h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasResonance+1, "N/A");
	}else{
		if (selection == "SingleMuon"){
			cout <<  "--- Setting up Single Muon Histograms --- " << endl;
			h_allMuonsPt = new TH1F("h_allMuonsPt", "h_allMuonsPt; Pt [GeV]; Muons", 250, 0,500);
			h_selectedMuonsPt = new TH1F("h_selectedMuonsPt", "h_selectedMuonsPt; Pt [GeV]; Muons", 250, 0,500);
			h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasResonance+1, "eventHasZ");
			float zInvMassMin = 0;
			float zInvMassMax = 130;
			h_osInvMass = new TH1F("h_osInvMass", "h_osInvMass; Mass [GeV]; Dimuons", invMassBins,zInvMassMin, zInvMassMax);
			h_ssInvMass = new TH1F("h_ssInvMass", "h_ssInvMass; Mass [GeV]; Dimuons", invMassBins, zInvMassMin, zInvMassMax);
			h_premassCutInvMass = new TH1F("h_premassCutInvMass", "h_preMassCutInvMass;Mass [GeV]; Dimuons", invMassBins,zInvMassMin, zInvMassMax);
		} else if (selection == "Charmonium") {
			cout <<  "--- Setting up Charmonium Histograms --- " << endl;
			h_allMuonsPt = new TH1F("h_allMuonsPt", "h_allMuonsPt; Pt [GeV]; Muons", 100, 0,200);
			h_selectedMuonsPt = new TH1F("h_selectedMuonsPt", "h_selectedMuonsPt; Pt [GeV]; Muons", 100, 0,200);
			h_eventCuts->GetXaxis()->SetBinLabel(EVENT_CUTS::hasResonance+1, "eventHasJ/#Psi");
			float jpsiMassMin = 2;
			float jpsiMassMax = 5;
			h_osInvMass = new TH1F("h_osInvMass", "h_osInvMass; Mass [GeV]; Dimuons", invMassBins,jpsiMassMin, jpsiMassMax);
			h_ssInvMass = new TH1F("h_ssInvMass", "h_ssInvMass; Mass [GeV]; Dimuons", invMassBins,jpsiMassMin, jpsiMassMax);
			h_premassCutInvMass = new TH1F("h_premassCutInvMass", "h_preMassCutInvMass;Mass [GeV]; Dimuons", invMassBins,jpsiMassMin, jpsiMassMax);
		}else if(selection == "MuonGun") {
			cout <<  "--- Tree instantiated as MuonGun sample --- " << endl;
			h_allMuonsPt = new TH1F("h_allMuonsPt", "h_allMuonsPt; Pt [GeV]; Muons", 400, 0,4000);
			h_selectedMuonsPt = new TH1F("h_selectedMuonsPt", "h_selectedMuonsPt; Pt [GeV]; Muons", 400, 0,4000);
		} else { //default to single muon selection
			cout <<  "--- Defaulting as singleMu sample --- " << endl;
		}
	}
}

SelectionHistograms::~SelectionHistograms(){
	if(h_eventCuts) delete h_eventCuts;
	if(h_muonCuts) delete h_muonCuts;
	if(h_osInvMass) delete h_osInvMass;
	if(h_ssInvMass) delete h_ssInvMass;
	if(h_premassCutInvMass) delete h_premassCutInvMass;
	if(h_nAllMuons) delete h_nAllMuons;
	if(h_allMuonsPt) delete h_allMuonsPt;
	if(h_allMuonsEta) delete h_allMuonsEta;
	if(h_allMuonsPhi) delete h_allMuonsPhi;
	if(h_allInvMass) delete h_allInvMass;
	if(h_selectedMuonsPt) delete h_selectedMuonsPt;
	if(h_selectedMuonsEta) delete h_selectedMuonsEta;
	if(h_selectedMuonsPhi) delete h_selectedMuonsPhi;
	if(h_nSelectedMuons) delete h_nSelectedMuons;
	if(h_nAllSegments) delete h_nAllSegments;
}

void SelectionHistograms::write(){
	if(h_eventCuts) h_eventCuts->Write();
	if(h_muonCuts) h_muonCuts->Write();
	if(h_osInvMass) h_osInvMass->Write();
	if(h_ssInvMass) h_ssInvMass->Write();
	if(h_premassCutInvMass) h_premassCutInvMass->Write();
	if(h_nAllMuons) h_nAllMuons->Write();
	if(h_allMuonsPt) h_allMuonsPt->Write();
	if(h_allMuonsEta) h_allMuonsEta->Write();
	if(h_allMuonsPhi) h_allMuonsPhi->Write();
	if(h_allInvMass) h_allInvMass->Write();
	if(h_selectedMuonsPt) h_selectedMuonsPt->Write();
	if(h_selectedMuonsEta) h_selectedMuonsEta->Write();
	if(h_selectedMuonsPhi) h_selectedMuonsPhi->Write();
	if(h_nSelectedMuons) h_nSelectedMuons->Write();
	if(h_nAllSegments) h_nAllSegments->Write();
}

void FillEventInfo::fill(const edm::Event& iEvent, unsigned int nSegments){
  EventNumber     	= iEvent.id().event();
  RunNumber       	= iEvent.id().run();
  LumiSection     	= iEvent.eventAuxiliary().luminosityBlock();
  BXCrossing      	= iEvent.eventAuxiliary().bunchCrossing();
  NSegmentsInEvent	= nSegments;
}

void FillGenParticleInfo::fill(const vector<reco::GenParticle>& gen){
	for(const auto& g :gen)
	{
		pdg_id	->push_back(g.pdgId());
		pt		->push_back(g.pt());
		eta		->push_back(g.eta());
		phi		->push_back(g.phi());
		q		->push_back(g.charge());
	}
}

void FillSimHitsInfo::fill(const vector<PSimHit>& simhits) {
	for(auto& sim : simhits)
	{
		CSCDetId id(sim.detUnitId());
		ch_id		->push_back(CSCHelper::serialize(id.station(),id.ring(), id.chamber(),id.endcap()));
		pdg_id		->push_back(sim.particleType());
		layer		->push_back(id.layer());
		energyLoss	->push_back(sim.energyLoss());
		thetaAtEntry->push_back(sim.thetaAtEntry());
		phiAtEntry	->push_back(sim.phiAtEntry());
		pAtEntry	->push_back(sim.pabs());
	}
}

void FillCaloHitsInfo::fill(const vector<PCaloHit>& calohits, const EcalBarrelGeometry* theEcal) 
{
	for(auto& cal : calohits) 
	{
		auto id = EBDetId(cal.id());
		if(!theEcal->present(id))
		{
			cout << "Barrel DetId not present" << endl;
		}

		auto geo = theEcal->getGeometry(id);
		if(!geo)
		{
			cout << "Not an ECAL!" << endl;
			return;
		}

		energyEM	->push_back(cal.energyEM());
		energyHad	->push_back(cal.energyHad());
		eta			->push_back(geo->etaPos());
		phi			->push_back(geo->phiPos());
	}
}

void FillCaloHitsInfo::fill(const vector<PCaloHit>& calohits, const EcalEndcapGeometry* theEcal) 
{
	for(auto& cal : calohits) 
	{
		auto id = EEDetId(cal.id());
		if(!theEcal->present(id))
		{
			cout << "Endcap DetId not present" << endl;
		}

		auto geo = theEcal->getGeometry(id);
		if(!geo){
			cout << "Not an ECAL!" << endl;
			return;
		}
		energyEM	->push_back(cal.energyEM());
		energyHad	->push_back(cal.energyHad());
		eta			->push_back(geo->etaPos());
		phi			->push_back(geo->phiPos());
	}
}

void FillCaloHitsInfo::fill(const vector<PCaloHit>& calohits, const EcalPreshowerGeometry* theEcal) {
	for(auto& cal : calohits) 
	{
		auto id = ESDetId(cal.id());
		if(!theEcal->present(id))
		{
			cout << "Preshower DetId not present" << endl;
		}

		auto geo = theEcal->getGeometry(id);
		if(!geo)
		{
			cout << "Not an ECAL!" << endl;
			return;
		}

		energyEM	->push_back(cal.energyEM());
		energyHad	->push_back(cal.energyHad());
		eta			->push_back(geo->etaPos());
		phi			->push_back(geo->phiPos());
	}
}

//TODO: could merge ECAL / HCAL fill functions
void FillCaloHitsInfo::fill(const vector<PCaloHit>& calohits, const HcalGeometry* theHcal)
{
	for(auto& cal: calohits){
		auto id = HcalDetId(cal.id());
		if(!theHcal->present(id)){
			cout << "HCAL DetId not present" << endl;
		}
		auto geo = theHcal->getGeometry(id);
		if(!geo) {
			cout << "Not an HCAL!" << endl;
			return;
		}
		energyEM	->push_back(cal.energyEM());
		energyHad	->push_back(cal.energyHad());
		eta			->push_back(geo->etaPos());
		phi			->push_back(geo->phiPos());
	}
}


void FillPFInfo::fill(const vector<reco::PFCandidate>& pfCand) 
{
	for(const auto& cand : pfCand) {
		pdg_id		->push_back(cand.pdgId());
		particleId	->push_back(cand.particleId());
		eta			->push_back(cand.eta());
		phi			->push_back(cand.phi());
		ecalEnergy	->push_back(cand.rawEcalEnergy());
		hcalEnergy	->push_back(cand.rawHcalEnergy());
		h0Energy	->push_back(cand.rawHoEnergy());
	}
}

void FillMuonInfo::fill(const reco::MuonCollection& muons)
{
	for (const auto &muon : muons) fill(muon);
}

void FillMuonInfo::fill(const reco::Muon& muon) {
	pt			->push_back(muon.pt());
	eta			->push_back(muon.eta());
	phi			->push_back(muon.phi());
	q			->push_back(muon.charge());
	isGlobal	->push_back(muon.isGlobalMuon());
	isTracker	->push_back(muon.isTrackerMuon());
}

size16 FillSegmentInfo::findRecHitIdx(const CSCRecHit2D& hit, const CSCRecHit2DCollection* allRecHits)
{
  int idx = -1;
  int foundIDX = -1;
  for (CSCRecHit2DCollection::const_iterator hiti=allRecHits->begin(); hiti!=allRecHits->end(); hiti++)
  {
    idx++;
    if (!hit.sharesInput(&(*hiti),CSCRecHit2D::all)) continue;
    foundIDX = idx;
    break;
  }
  if(foundIDX < 0) throw std::invalid_argument("FillSegmentInfo::findRecHitIdx -> Could not find rechit");
  return CSCHelper::convertTo<size16,int>(foundIDX,"foundIDX");
}


void FillSegmentInfo::fill(std::vector<const CSCSegment*>& segments, const CSCGeometry* theCSC, int mu_index){
  for (auto dSiter=segments.begin(); dSiter != segments.end(); dSiter++) 
  {
	 fill(**dSiter, theCSC, mu_index);
  }
}

void FillSegmentInfo::fill(const CSCSegment& segment, const CSCGeometry* theCSC, int mu_index)
{
	DetId detId  = segment.geographicalId();
	CSCDetId id(detId.rawId());

	const auto& segmentHits = segment.specificRecHits();
	//const auto& segmentHits = dSiter->specificRecHits();
	//inherited from C. Bravo - translates [cm] -> [strips]
	const CSCChamber *segChamber = theCSC->chamber(id);
	const CSCLayer *segLay3 = segChamber->layer(3);
	const CSCLayer *segLay4 = segChamber->layer(4);
	const CSCLayerGeometry *segLay3Geo = segLay3->geometry();
	LocalPoint segLPlayer = segLay3->toLocal(segChamber->toGlobal(segment.localPosition()));
	LocalVector segLVlayer = segLay3->toLocal(segChamber->toGlobal(segment.localDirection()));
	float scale = -1.0*segLPlayer.z()/segLVlayer.z();
	LocalVector tV = scale*segLVlayer;
	LocalPoint tP = segLPlayer + tV;
	float segStr = segLay3Geo->strip(tP);
	int strI = floor(segStr);
	float cm2strip = fabs( segLay3Geo->xOfStrip(strI,tP.y()) - segLay3Geo->xOfStrip(strI+1,tP.y()) );

	LocalPoint lzero(0.0,0.0,0.0);
	GlobalPoint lay4zero = segLay4->toGlobal(lzero);
	LocalPoint lay4zeroIn3 = segLay3->toLocal(lay4zero);
	float cm2lay = fabs(lay4zeroIn3.z());


	//fill everything
	mu_id	->push_back(mu_index);
	ch_id	->push_back(CSCHelper::serialize(id.station(), id.ring(), id.chamber(), id.endcap()));
	pos_x	->push_back(segLay3Geo->strip(tP));
	pos_y	->push_back(segLay3Geo->wireGroup(segLay3Geo->nearestWire(tP)));
	dxdz	->push_back(segment.localDirection().x() / cm2strip / ( segment.localDirection().z() / cm2lay));
	dydz	->push_back(segment.localDirection().y() / segment.localDirection().z());
	chisq	->push_back(segment.chi2());
	nHits	->push_back(CSCHelper::convertTo<size8>(segmentHits.size(), "segment_nHits"));

}

void FillRecHitInfo::fill(const std::vector<CSCRecHit2D>& recHits, int mu_index)
{
  for (const auto& rechit : recHits){

      DetId detId = rechit.geographicalId();
      CSCDetId id(detId.rawId());

      int centerID = rechit.nStrips()/2;
      int centerStr = rechit.channels(centerID);

      float rhMaxBuf = -999.0;
      for(int tI = 0; tI < int(rechit.nTimeBins()); tI++)
      {
    	  if(rechit.adcs(centerID,tI) > rhMaxBuf) rhMaxBuf = rechit.adcs(centerID,tI);
      }

      mu_id->push_back(mu_index);
      ch_id->push_back(CSCHelper::serialize(id.station(), id.ring(), id.chamber(), id.endcap()));
      lay->push_back(CSCHelper::convertTo<size8>(id.layer(),"rh_lay"));
      pos_x->push_back(float(centerStr) + rechit.positionWithinStrip());
      pos_y->push_back(-1); //NOT IMPLEMENTED
      e->push_back(rechit.energyDepositedInLayer());
      max_adc->push_back(rhMaxBuf);
  }
}

void FillLCTInfo::fill(const CSCCorrelatedLCTDigiCollection& lcts)
{
  for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator chamber=lcts.begin(); chamber!=lcts.end(); chamber++)
  {
    CSCDetId id = (*chamber).first;
    unsigned int st = id.station();
    unsigned int ri = id.ring();
    unsigned int ch = id.chamber();
    unsigned int ec = id.endcap();

    const CSCCorrelatedLCTDigiCollection::Range& range =(*chamber).second;
    for(CSCCorrelatedLCTDigiCollection::const_iterator digiItr = range.first; digiItr != range.second; ++digiItr)
    {

    	if(st == 1 && (ri == 1|| ri ==4)){
    		//we need to manually adjust this because they don't for us
    		//getStrip returns a half strip (who knows why)
    		if(digiItr->getStrip() > CSCHelper::MAX_ME11B_HALF_STRIP) ri = 4;
    		else ri = 1; //resets ring in case where multiple clcts in ME11
    	}

      ch_id       	->push_back(CSCHelper::serialize(st, ri, ch, ec));
      quality     	->push_back(CSCHelper::convertTo<size8>(digiItr->getQuality(),"lct_quality"));
      pattern     	->push_back(CSCHelper::convertTo<size8>(digiItr->getPattern(),"lct_pattern"));
      bend        	->push_back(CSCHelper::convertTo<size8>(digiItr->getBend()   ,"lct_bend"));
      keyWireGroup	->push_back(CSCHelper::convertTo<size8>(digiItr->getKeyWG()  ,"lct_keyWireGroup"));
      keyHalfStrip	->push_back(CSCHelper::convertTo<size8>(digiItr->getStrip()  ,"lct_keyHalfStrip"));
	  bunchCross	->push_back(CSCHelper::convertTo<size8>(digiItr->getBX()	 ,"lct_bunchCross"));
      }
  }
}

void FillCLCTInfo::fill(const CSCCLCTDigiCollection& clcts) 
{
	for (CSCCLCTDigiCollection::DigiRangeIterator chamber = clcts.begin();
			chamber != clcts.end(); chamber++) {
		CSCDetId id = (*chamber).first;
	    unsigned int st = id.station();
	    unsigned int ri = id.ring();
	    unsigned int ch = id.chamber();
	    unsigned int ec = id.endcap();

		const CSCCLCTDigiCollection::Range& range = (*chamber).second;
		for (CSCCLCTDigiCollection::const_iterator digiItr = range.first;
				digiItr != range.second; ++digiItr) {

			if(st == 1 && (ri == 1|| ri ==4)){
				//we need to manually adjust this because they don't for us
				//getStrip returns a half strip (who knows why)
				if(digiItr->getKeyStrip() > CSCHelper::MAX_ME11B_HALF_STRIP) ri = 4;
				else ri = 1; //resets ring in case where multiple clcts in ME11
			}
			ch_id->push_back(CSCHelper::serialize(st, ri, ch, ec));
			isValid->push_back(
				CSCHelper::convertTo<size8>(digiItr->isValid(),
						"clct_isvalid"));
			quality->push_back(
				CSCHelper::convertTo<size8>(digiItr->getQuality(),
						"clct_quality"));
			pattern->push_back(
				CSCHelper::convertTo<size8>(digiItr->getPattern(),
						"clct_pattern"));
			stripType->push_back(
				CSCHelper::convertTo<size8>(digiItr->getStripType(),
						"clct_stripType"));
			bend->push_back(
				CSCHelper::convertTo<size8>(digiItr->getBend(),
						"clct_bend"));
			halfStrip->push_back(
				CSCHelper::convertTo<size8>(digiItr->getStrip(),
						"clct_halfStrip"));
			CFEB->push_back(
				CSCHelper::convertTo<size8>(digiItr->getCFEB(),
						"clct_CFEB"));
			BX->push_back(
				CSCHelper::convertTo<size8>(digiItr->getBX(), "clct_BX"));
			trkNumber->push_back(
				CSCHelper::convertTo<size8>(digiItr->getTrknmb(),
						"clct_trkNumber"));
			keyStrip->push_back(
					CSCHelper::convertTo<size8>(digiItr->getKeyStrip(),
							"clct_keyStrip"));
		}

	}
}

void FillALCTInfo::fill(const CSCALCTDigiCollection &alcts)
{
	CSCALCTDigiCollection::DigiRangeIterator chamber;
	for (chamber = alcts.begin(); chamber!=alcts.end(); chamber++)
	{
		CSCDetId id = (*chamber).first;
		unsigned int st = id.station();
		unsigned int ri = id.ring();
		unsigned int ch = id.chamber();
	  	unsigned int ec = id.endcap();

		const CSCALCTDigiCollection::Range &range = (*chamber).second;
		CSCALCTDigiCollection::const_iterator digiItr;

		for (digiItr = range.first; digiItr!=range.second; digiItr++)
		{
			ch_id		->push_back(CSCHelper::serialize(st, ri, ch, ec));
			isValid		->push_back(CSCHelper::convertTo<size8>(digiItr->isValid(),"alct_isValid"));
			quality		->push_back(CSCHelper::convertTo<size8>(digiItr->getQuality(),"alct_quality"));
			accelerator	->push_back(CSCHelper::convertTo<size8>(digiItr->getAccelerator(),"alct_accelerator"));
			collisionB	->push_back(CSCHelper::convertTo<size8>(digiItr->getCollisionB(),"alct_collisionB"));
			keyWG		->push_back(CSCHelper::convertTo<size8>(digiItr->getKeyWG(),"alct_keyWG"));
			BX			->push_back(CSCHelper::convertTo<size8>(digiItr->getBX(),"alct_BX"));
			trkNumber	->push_back(CSCHelper::convertTo<size8>(digiItr->getTrknmb(),"alct_trkNumber"));
			fullBX		->push_back(CSCHelper::convertTo<size16>(digiItr->getFullBX(),"alct_fullBX"));
		}
	}
}

void FillWireInfo::fill(const CSCWireDigiCollection &wires)
{
	CSCWireDigiCollection::DigiRangeIterator chamber;
	for (chamber = wires.begin(); chamber!= wires.end(); chamber++)
	{
		CSCDetId id = (*chamber).first;
		unsigned int st = id.station();
		unsigned int ri = id.ring();
		unsigned int ch = id.chamber();
		unsigned int ec = id.endcap();

		const CSCWireDigiCollection::Range &range = (*chamber).second;
		CSCWireDigiCollection::const_iterator digiItr;

		for (digiItr = range.first; digiItr!= range.second; digiItr++)
		{
			ch_id			->push_back(CSCHelper::serialize(st, ri, ch, ec));
			group			->push_back(digiItr->getWireGroup());					
			lay				->push_back(CSCHelper::convertTo<size8>(id.layer(), "wire_lay"));
			timeBin			->push_back(digiItr->getTimeBin());
			timeBinWord		->push_back(digiItr->getTimeBinWord());
			BX				->push_back(digiItr->getWireGroupBX());
		}
	}
}

void FillStripInfo::fill(const CSCStripDigiCollection &strips)
{
	CSCStripDigiCollection::DigiRangeIterator chamber;
	for (chamber = strips.begin(); chamber != strips.end(); chamber++)
	{
		CSCDetId id = (*chamber).first;
		unsigned int st = id.station();
		unsigned int ri = id.ring();
		unsigned int ch = id.chamber();
		unsigned int ec = id.endcap();

		const CSCStripDigiCollection::Range &range = (*chamber).second;
		CSCStripDigiCollection::const_iterator digiItr;

		for (digiItr = range.first; digiItr != range.second; digiItr++)
		{
			if (st == 1 && (ri == 1|| ri ==4))
			{
				// we need to manually adjust this because they don't for us
				// getStrip returns a strip this time (different than before)

				if(digiItr->getStrip() > CSCHelper::MAX_ME11B_STRIP) ri = 4;
				else ri = 1; // resets ring in case where multiple clcts in ME11
			}

			ch_id				->push_back(CSCHelper::serialize(st, ri, ch, ec));
			lay					->push_back(CSCHelper::convertTo<size8>(id.layer(), "strip_lay"));
			num					->push_back(digiItr->getStrip());
			//ADC				->push_back(digiItr->getADCCounts());
			//L1APhase			->push_back(digiItr->getL1APhase())
			//ADCOverflow		->push_back(digiItr->getADCOverflow());
			//OverlappedSample	->push_back(digiItr->getOverlappedSample());
			//Errorstat			->push_back(digiItr->getErrorstat());
		}
	}
}

void FillCompInfo::fill(const CSCComparatorDigiCollection& comps){
  for (CSCComparatorDigiCollection::DigiRangeIterator chamber=comps.begin(); chamber!=comps.end(); chamber++)
  {
    CSCDetId id = (*chamber).first;
    unsigned int st = id.station();
    unsigned int ri = id.ring();
    unsigned int ch = id.chamber();
    unsigned int ec = id.endcap();
    const CSCComparatorDigiCollection::Range& range =(*chamber).second;
    for(CSCComparatorDigiCollection::const_iterator digiItr = range.first; digiItr != range.second; ++digiItr)
    {

    	if(st == 1 && (ri == 1|| ri ==4))
		{
    		//we need to manually adjust this because they don't for us
    		//getStrip returns a strip this time (different than before)
			
    		if(digiItr->getStrip() > CSCHelper::MAX_ME11B_STRIP) ri = 4;
    		else ri = 1; //resets ring in case where multiple clcts in ME11
    	}
    	ch_id       ->push_back(CSCHelper::serialize(st, ri, ch, ec));
		lay			->push_back(CSCHelper::convertTo<size8>(id.layer(), "comp_lay"));
		strip		->push_back(CSCHelper::convertTo<size8>((*digiItr).getStrip(),"comp_strip"));
		halfStrip	->push_back(CSCHelper::convertTo<size8>((*digiItr).getComparator(),"comp_comp"));
		nTimeOn		->push_back(CSCHelper::convertTo<size8>((*digiItr).getTimeBinsOn().size(),"comp_nTimeOn"));
		unsigned int bestTimeBin = 0;
		for(const auto& time : (*digiItr).getTimeBinsOn())
		{
			bestTimeBin = time;
			if (time < 10 && time > 5) break;
		}
		bestTime	->push_back(CSCHelper::convertTo<size8>(bestTimeBin, "comp_bestTime"));
		}
	}
}

/*
void FillStripInfo::fill(const CSCStripDigiCollection& strips){
  reset();
  for (CSCStripDigiCollection::DigiRangeIterator dSDiter=strips.begin(); dSDiter!=strips.end(); dSDiter++) {
    CSCDetId id = (CSCDetId)(*dSDiter).first;

    std::vector<CSCStripDigi>::const_iterator stripIter = (*dSDiter).second.first;
    std::vector<CSCStripDigi>::const_iterator lStrip = (*dSDiter).second.second;
    for( ; stripIter != lStrip; ++stripIter) {
      std::vector<int> myADCVals = stripIter->getADCCounts();
      bool thisStripFired = false;
      float thisPedestal = 0.5*(float)(myADCVals[0]+myADCVals[1]);
      float threshold = 13.3 ;
      float diff = 0.;
      for (unsigned int iCount = 0; iCount < myADCVals.size(); iCount++) {
        diff = (float)myADCVals[iCount]-thisPedestal;
        if (diff > threshold) { thisStripFired = true; }
      }
      if(!thisStripFired) continue;

      strip_id.push_back(CSCHelper::chamberSerial(id));
      strip_lay.push_back(CSCHelper::convertTo<size8>(id.layer(),"strip_lay"));
      strip_number.push_back(CSCHelper::convertTo<size8>(stripIter->getStrip(),"strip_number"));
    }
  } // end strip loop
}
*/










