#ifndef Phase2OuterTracker_OuterTrackerTrack_h
#define Phase2OuterTracker_OuterTrackerTrack_h

#include <vector>
#include <memory>
//#include "FWCore/Utilities/interface/EDGetToken.h"
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
//#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"  // does not exist yet in 90X
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"  // does not exist yet in 90X

class DQMStore;

class OuterTrackerTrack : public edm::EDAnalyzer {

public:
  explicit OuterTrackerTrack(const edm::ParameterSet&);
  ~OuterTrackerTrack();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  
  
  MonitorElement* Track_NStubs = 0;
  
  /// Low-quality TTTracks (made from less than X TTStubs)
  MonitorElement* Track_LQ_N = 0;
  MonitorElement* Track_LQ_Pt = 0;
  MonitorElement* Track_LQ_Eta = 0;
  MonitorElement* Track_LQ_Phi = 0;
  MonitorElement* Track_LQ_VtxZ0 = 0;
  MonitorElement* Track_LQ_Chi2 = 0;
  MonitorElement* Track_LQ_Chi2Red = 0;
  MonitorElement* Track_LQ_Chi2_NStubs = 0;
  MonitorElement* Track_LQ_Chi2Red_NStubs = 0;
  
  /// High-quality TTTracks (made from at least X TTStubs)
  MonitorElement* Track_HQ_N = 0;
  MonitorElement* Track_HQ_Pt = 0;
  MonitorElement* Track_HQ_Eta = 0;
  MonitorElement* Track_HQ_Phi = 0;
  MonitorElement* Track_HQ_VtxZ0 = 0;
  MonitorElement* Track_HQ_Chi2 = 0;
  MonitorElement* Track_HQ_Chi2Red = 0;
  MonitorElement* Track_HQ_Chi2_NStubs = 0;
  MonitorElement* Track_HQ_Chi2Red_NStubs = 0;
  
  
 private:
  DQMStore* dqmStore_;
  edm::ParameterSet conf_;
  //edm::EDGetTokenT<edmNew::DetSetVector< TTTrack< Ref_Phase2TrackerDigi_ > > >  tagTTTracksToken_;
  //edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >  tagTTTrackMCTruthToken_;

  std::string topFolderName_;
  unsigned int HQDelim_;
  bool verbosePlots_;
  
};
#endif
