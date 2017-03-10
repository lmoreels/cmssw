#ifndef Phase2OuterTracker_OuterTrackerCluster_h
#define Phase2OuterTracker_OuterTrackerCluster_h

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
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"

class DQMStore;

class OuterTrackerCluster : public edm::EDAnalyzer {

public:
  explicit OuterTrackerCluster(const edm::ParameterSet&);
  ~OuterTrackerCluster();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //virtual void beginJob() ;
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);
  
  // TTCluster stacks
  MonitorElement* Cluster_Gen_Barrel = 0;
  MonitorElement* Cluster_Unkn_Barrel = 0;
  MonitorElement* Cluster_Comb_Barrel = 0;
  MonitorElement* Cluster_Gen_Endcap_Disc = 0;
  MonitorElement* Cluster_Unkn_Endcap_Disc = 0;
  MonitorElement* Cluster_Comb_Endcap_Disc = 0;
  MonitorElement* Cluster_Gen_Endcap_Ring = 0;
  MonitorElement* Cluster_Unkn_Endcap_Ring = 0;
  MonitorElement* Cluster_Comb_Endcap_Ring = 0;
  MonitorElement* Cluster_Gen_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Unkn_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Comb_Endcap_Ring_Fw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Gen_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Unkn_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Comb_Endcap_Ring_Bw[5] = {0, 0, 0, 0, 0};
  MonitorElement* Cluster_Gen_Eta = 0;
  MonitorElement* Cluster_Unkn_Eta = 0;
  MonitorElement* Cluster_Comb_Eta = 0;
  
  
 private:
  DQMStore* dqmStore_;
  edm::ParameterSet conf_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > >  tagTTClustersToken_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > >  tagTTClusterMCTruthToken_;

  std::string topFolderName_;
  int nDiscs_;
  bool verbosePlots_;
  
};
#endif
