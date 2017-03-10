#ifndef Phase2OuterTracker_OuterTrackerStub_h
#define Phase2OuterTracker_OuterTrackerStub_h

#include <vector>
#include <memory>
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DQM/SiStripCommon/interface/TkHistoMap.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"	//Needed for TTStubAssociationMap.h !
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"


class DQMStore;

class OuterTrackerStub : public edm::EDAnalyzer {

public:
  explicit OuterTrackerStub(const edm::ParameterSet&);
  ~OuterTrackerStub();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob() ;
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  
  // TTStub stacks
  MonitorElement* Stub_Gen_Barrel = 0;
  MonitorElement* Stub_Unkn_Barrel = 0;
  MonitorElement* Stub_Comb_Barrel = 0;
  MonitorElement* Stub_Gen_Endcap_Disc = 0;
  MonitorElement* Stub_Unkn_Endcap_Disc = 0;
  MonitorElement* Stub_Comb_Endcap_Disc = 0;
  MonitorElement* Stub_Gen_Endcap_Ring = 0;
  MonitorElement* Stub_Unkn_Endcap_Ring = 0;
  MonitorElement* Stub_Comb_Endcap_Ring = 0;
  MonitorElement* Stub_Gen_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Unkn_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Comb_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Gen_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Unkn_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Comb_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Stub_Gen_Eta = 0;
  MonitorElement* Stub_Unkn_Eta = 0;
  MonitorElement* Stub_Comb_Eta = 0;
  
  
 private:
  DQMStore* dqmStore_;
  edm::ParameterSet conf_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >  tagTTStubsToken_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > >  tagTTStubMCTruthToken_;

  std::string topFolderName_;
  int nDiscs_;
  bool verbosePlots_;
  
};
#endif
