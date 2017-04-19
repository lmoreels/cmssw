// -*- C++ -*-
//
// Package:    Phase2OuterTracker
// Class:      Phase2OuterTracker
// 
/**\class Phase2OuterTracker OuterTrackerCluster.cc Validation/Phase2OuterTracker/plugins/OuterTrackerCluster.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lieselotte Moreels
//         Created:  Fri, 13 Jun 2014 09:57:34 GMT
// $Id$
//
//

// system include files
#include <memory>
#include <vector>
#include <numeric>
#include <fstream>
#include <math.h>
#include "TNamed.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Validation/Phase2OuterTracker/interface/OuterTrackerCluster.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "TMath.h"
#include <iostream>

//
// constructors and destructor
//
OuterTrackerCluster::OuterTrackerCluster(const edm::ParameterSet& iConfig)
: dqmStore_(edm::Service<DQMStore>().operator->()), conf_(iConfig)
{
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  tagTTClustersToken_ = consumes<edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > (conf_.getParameter<edm::InputTag>("TTClusters") );
  tagTTClusterMCTruthToken_ = consumes< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > (conf_.getParameter<edm::InputTag>("TTClusterMCTruth") );
  nDiscs_ = conf_.getUntrackedParameter<int>("nDiscs", 5);
  verbosePlots_ = conf_.getUntrackedParameter<bool>("verbosePlots",false);
}


OuterTrackerCluster::~OuterTrackerCluster()
{
	
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
	
}


//
// member functions
//

// ------------ method called for each event  ------------
void
OuterTrackerCluster::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  /// Track Trigger
  edm::Handle< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > > Phase2TrackerDigiTTClusterHandle;
  iEvent.getByToken( tagTTClustersToken_, Phase2TrackerDigiTTClusterHandle );
  /// Track Trigger MC Truth
  edm::Handle< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTClusterHandle;
  iEvent.getByToken( tagTTClusterMCTruthToken_, MCTruthTTClusterHandle );
  
  if ( ! Phase2TrackerDigiTTClusterHandle.isValid() || ! MCTruthTTClusterHandle.isValid() )  return;
  
  /// Geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  const TrackerTopology* tTopo;
  iSetup.get< TrackerTopologyRcd >().get(tTopoHandle);
  tTopo = tTopoHandle.product();
  
  edm::ESHandle< TrackerGeometry > tGeometryHandle;
  const TrackerGeometry* theTrackerGeometry;
  iSetup.get< TrackerDigiGeometryRecord >().get( tGeometryHandle );
  theTrackerGeometry = tGeometryHandle.product();
  
  	
  /// Loop over the input Clusters
  typename edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  typename edmNew::DetSet< TTCluster< Ref_Phase2TrackerDigi_ > >::const_iterator contentIter;
  for ( inputIter = Phase2TrackerDigiTTClusterHandle->begin();
       inputIter != Phase2TrackerDigiTTClusterHandle->end();
       ++inputIter )
  {
    for ( contentIter = inputIter->begin();
         contentIter != inputIter->end();
         ++contentIter )
    {
      /// Make the reference to be put in the map
      edm::Ref< edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > >, TTCluster< Ref_Phase2TrackerDigi_ > > tempCluRef = edmNew::makeRefTo( Phase2TrackerDigiTTClusterHandle, contentIter );
      
      DetId detIdClu = theTrackerGeometry->idToDet( tempCluRef->getDetId() )->geographicalId();
      
      bool genuineClu     = MCTruthTTClusterHandle->isGenuine( tempCluRef );
      bool combinClu      = MCTruthTTClusterHandle->isCombinatoric( tempCluRef );
      
      
      MeasurementPoint mp = tempCluRef->findAverageLocalCoordinates();
      const GeomDet* theGeomDet = theTrackerGeometry->idToDet(detIdClu);
      Global3DPoint posClu = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(mp) );
      
      if ( detIdClu.subdetId() == static_cast<int>(StripSubdetector::TOB) )  // Phase 2 Outer Tracker Barrel
      {
        if ( genuineClu )
        {
          Cluster_Gen_Barrel->Fill( tTopo->layer(detIdClu) );
        }
        else if ( combinClu )
        {
          Cluster_Comb_Barrel->Fill( tTopo->layer(detIdClu) );
        }
        else
        {
          Cluster_Unkn_Barrel->Fill( tTopo->layer(detIdClu) );
        }
      }	// end if isBarrel
      else if ( detIdClu.subdetId() == static_cast<int>(StripSubdetector::TID) )  // Phase 2 Outer Tracker Endcap
      {
        int side = tTopo->side(detIdClu);
        int disc = tTopo->layer(detIdClu);  // returns wheel
        int ring = tTopo->tidRing(detIdClu);
        if ( genuineClu )
        {
          Cluster_Gen_Endcap_Disc->Fill( disc );
          Cluster_Gen_Endcap_Ring->Fill( ring );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Cluster_Gen_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Cluster_Gen_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
        else if ( combinClu )
        {
          Cluster_Comb_Endcap_Disc->Fill( disc );
          Cluster_Comb_Endcap_Ring->Fill( ring );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Cluster_Comb_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Cluster_Comb_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
        else
        {
          Cluster_Unkn_Endcap_Disc->Fill( disc );
          Cluster_Unkn_Endcap_Ring->Fill( tTopo->tidRing(detIdClu) );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Cluster_Unkn_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Cluster_Unkn_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
      }	// end if isEndcap
      
      /// Eta distribution in function of genuine/combinatorial/unknown cluster
      if ( genuineClu ) Cluster_Gen_Eta->Fill( posClu.eta() );
      else if ( combinClu ) Cluster_Comb_Eta->Fill( posClu.eta() );
      else Cluster_Unkn_Eta->Fill( posClu.eta() );
    } // end loop contentIter
  } // end loop inputIter
  
}


// ------------ method called once each job just before starting event loop  ------------
void
OuterTrackerCluster::beginRun(const edm::Run& run, const edm::EventSetup& es)
{
  std::string HistoName;
  
  dqmStore_->setCurrentFolder(topFolderName_+"/TTClusters/");
  
  edm::ParameterSet psTTClusterEta =  conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Eta");
  HistoName = "Cluster_Gen_Eta";
  Cluster_Gen_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterEta.getParameter<int32_t>("Nbinsx"),
      psTTClusterEta.getParameter<double>("xmin"),
      psTTClusterEta.getParameter<double>("xmax"));
  Cluster_Gen_Eta->setAxisTitle("#eta", 1);
  Cluster_Gen_Eta->setAxisTitle("# Genuine L1 Clusters", 2);
  
  HistoName = "Cluster_Unkn_Eta";
  Cluster_Unkn_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterEta.getParameter<int32_t>("Nbinsx"),
      psTTClusterEta.getParameter<double>("xmin"),
      psTTClusterEta.getParameter<double>("xmax"));
  Cluster_Unkn_Eta->setAxisTitle("#eta", 1);
  Cluster_Unkn_Eta->setAxisTitle("# Unknown L1 Clusters", 2);
  
  HistoName = "Cluster_Comb_Eta";
  Cluster_Comb_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterEta.getParameter<int32_t>("Nbinsx"),
      psTTClusterEta.getParameter<double>("xmin"),
      psTTClusterEta.getParameter<double>("xmax"));
  Cluster_Comb_Eta->setAxisTitle("#eta", 1);
  Cluster_Comb_Eta->setAxisTitle("# Combinatorial L1 Clusters", 2);
  
  
  /// TTCluster stacks
  edm::ParameterSet psTTClusterLayers =  conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Layers");
  HistoName = "NClusters_Gen_Barrel";
  Cluster_Gen_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterLayers.getParameter<int32_t>("Nbinsx"),
      psTTClusterLayers.getParameter<double>("xmin"),
      psTTClusterLayers.getParameter<double>("xmax"));
  Cluster_Gen_Barrel->setAxisTitle("Barrel Layer", 1);
  Cluster_Gen_Barrel->setAxisTitle("# Genuine L1 Clusters", 2);
  
  HistoName = "NClusters_Unkn_Barrel";
  Cluster_Unkn_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterLayers.getParameter<int32_t>("Nbinsx"),
      psTTClusterLayers.getParameter<double>("xmin"),
      psTTClusterLayers.getParameter<double>("xmax"));
  Cluster_Unkn_Barrel->setAxisTitle("Barrel Layer", 1);
  Cluster_Unkn_Barrel->setAxisTitle("# Unknown L1 Clusters", 2);
  
  HistoName = "NClusters_Comb_Barrel";
  Cluster_Comb_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterLayers.getParameter<int32_t>("Nbinsx"),
      psTTClusterLayers.getParameter<double>("xmin"),
      psTTClusterLayers.getParameter<double>("xmax"));
  Cluster_Comb_Barrel->setAxisTitle("Barrel Layer", 1);
  Cluster_Comb_Barrel->setAxisTitle("# Combinatorial L1 Clusters", 2);
  
  edm::ParameterSet psTTClusterDiscs =  conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Discs");
  HistoName = "NClusters_Gen_Endcap_Disc";
  Cluster_Gen_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterDiscs.getParameter<int32_t>("Nbinsx"),
      psTTClusterDiscs.getParameter<double>("xmin"),
      psTTClusterDiscs.getParameter<double>("xmax"));
  Cluster_Gen_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Cluster_Gen_Endcap_Disc->setAxisTitle("# Genuine L1 Clusters", 2);
  
  HistoName = "NClusters_Unkn_Endcap_Disc";
  Cluster_Unkn_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterDiscs.getParameter<int32_t>("Nbinsx"),
      psTTClusterDiscs.getParameter<double>("xmin"),
      psTTClusterDiscs.getParameter<double>("xmax"));
  Cluster_Unkn_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Cluster_Unkn_Endcap_Disc->setAxisTitle("# Unknown L1 Clusters", 2);
  
  HistoName = "NClusters_Comb_Endcap_Disc";
  Cluster_Comb_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterDiscs.getParameter<int32_t>("Nbinsx"),
      psTTClusterDiscs.getParameter<double>("xmin"),
      psTTClusterDiscs.getParameter<double>("xmax"));
  Cluster_Comb_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Cluster_Comb_Endcap_Disc->setAxisTitle("# Combinatorial L1 Clusters", 2);
  
  edm::ParameterSet psTTClusterRings =  conf_.getParameter<edm::ParameterSet>("TH1TTCluster_Rings");
  HistoName = "NClusters_Gen_Endcap_Ring";
  Cluster_Gen_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterRings.getParameter<int32_t>("Nbinsx"),
      psTTClusterRings.getParameter<double>("xmin"),
      psTTClusterRings.getParameter<double>("xmax"));
  Cluster_Gen_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Cluster_Gen_Endcap_Ring->setAxisTitle("# Genuine L1 Clusters", 2);
  
  HistoName = "NClusters_Unkn_Endcap_Ring";
  Cluster_Unkn_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterRings.getParameter<int32_t>("Nbinsx"),
      psTTClusterRings.getParameter<double>("xmin"),
      psTTClusterRings.getParameter<double>("xmax"));
  Cluster_Unkn_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Cluster_Unkn_Endcap_Ring->setAxisTitle("# Unknown L1 Clusters", 2);
  
  HistoName = "NClusters_Comb_Endcap_Ring";
  Cluster_Comb_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTClusterRings.getParameter<int32_t>("Nbinsx"),
      psTTClusterRings.getParameter<double>("xmin"),
      psTTClusterRings.getParameter<double>("xmax"));
  Cluster_Comb_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Cluster_Comb_Endcap_Ring->setAxisTitle("# Combinatorial L1 Clusters", 2);
  
  /// Plots for debugging
  if ( verbosePlots_ )
  {
    dqmStore_->setCurrentFolder(topFolderName_+"/TTClusters/NClustersPerRing");
    
    for (int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Gen_Disc+"+std::to_string(i+1);
      Cluster_Gen_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Gen_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Gen_Endcap_Ring_Fw[i]->setAxisTitle("# Genuine L1 Clusters", 2);
    }

    for (int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Gen_Disc-"+std::to_string(i+1);
      Cluster_Gen_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Gen_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Gen_Endcap_Ring_Bw[i]->setAxisTitle("# Genuine L1 Clusters", 2);
    }

    for (int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Unkn_Disc+"+std::to_string(i+1);
      Cluster_Unkn_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("# Unknown L1 Clusters", 2);
    }

    for(int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Unkn_Disc-"+std::to_string(i+1);
      Cluster_Unkn_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("# Unknown L1 Clusters", 2);
    }

    for (int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Comb_Disc+"+std::to_string(i+1);
      Cluster_Comb_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Comb_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Comb_Endcap_Ring_Fw[i]->setAxisTitle("# Combinatorial L1 Clusters", 2);
    }

    for (int i = 0; i < nDiscs_; i++)
    {
      HistoName = "NClusters_Comb_Disc-"+std::to_string(i+1);
      Cluster_Comb_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTClusterRings.getParameter<int32_t>("Nbinsx"),
          psTTClusterRings.getParameter<double>("xmin"),
          psTTClusterRings.getParameter<double>("xmax"));
      Cluster_Comb_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Cluster_Comb_Endcap_Ring_Bw[i]->setAxisTitle("# Combinatorial L1 Clusters", 2);
    }
  
  }  /// end verbosePlots
  
} //end of method

// ------------ method called once each job just after ending the event loop  ------------
void 
OuterTrackerCluster::endJob(void) 
{
	
}

//define this as a plug-in
DEFINE_FWK_MODULE(OuterTrackerCluster);
