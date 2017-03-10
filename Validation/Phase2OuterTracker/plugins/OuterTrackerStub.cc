// -*- C++ -*-
//
// Package:    Phase2OuterTracker
// Class:      Phase2OuterTracker
// 
/**\class Phase2OuterTracker OuterTrackerStub.cc Validation/Phase2OuterTracker/plugins/OuterTrackerStub.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lieselotte Moreels
//         Created:  Mon, 27 Oct 2014 09:07:51 GMT
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
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "Validation/Phase2OuterTracker/interface/OuterTrackerStub.h"
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
OuterTrackerStub::OuterTrackerStub(const edm::ParameterSet& iConfig)
: dqmStore_(edm::Service<DQMStore>().operator->()), conf_(iConfig)

{
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  tagTTStubsToken_ = consumes<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > (conf_.getParameter<edm::InputTag>("TTStubs") );
  tagTTStubMCTruthToken_ = consumes<edmNew::DetSetVector< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > > (conf_.getParameter<edm::InputTag>("TTStubMCTruth") );
  verbosePlots_ = conf_.getUntrackedParameter<bool>("verbosePlots",false);
}

OuterTrackerStub::~OuterTrackerStub()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}

//
// member functions
//

// ------------ method called for each event  ------------
void
OuterTrackerStub::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /// Track Trigger
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > Phase2TrackerDigiTTStubHandle;
  iEvent.getByToken( tagTTStubsToken_, Phase2TrackerDigiTTStubHandle );
  /// Track Trigger MC Truth
  edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
  iEvent.getByToken( tagTTStubMCTruthToken_, MCTruthTTStubHandle );
  
  /// Geometry
  edm::ESHandle<TrackerTopology> tTopoHandle;
  const TrackerTopology* tTopo;
  iSetup.get< TrackerTopologyRcd >().get(tTopoHandle);
  tTopo = tTopoHandle.product();
  
  edm::ESHandle< TrackerGeometry > tGeometryHandle;
  const TrackerGeometry* theTrackerGeometry;
  iSetup.get< TrackerDigiGeometryRecord >().get( tGeometryHandle );
  theTrackerGeometry = tGeometryHandle.product();
  
  
  /// Loop over the input Stubs
  typename edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator inputIter;
  typename edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > >::const_iterator contentIter;
  for ( inputIter = Phase2TrackerDigiTTStubHandle->begin();
       inputIter != Phase2TrackerDigiTTStubHandle->end();
       ++inputIter )
  {
    for ( contentIter = inputIter->begin();
         contentIter != inputIter->end();
         ++contentIter )
    {
      /// Make the reference to be put in the map
      edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > >, TTStub< Ref_Phase2TrackerDigi_ > > tempStubRef = edmNew::makeRefTo( Phase2TrackerDigiTTStubHandle, contentIter );
      
      /// Get det ID (place of the stub)
      //  tempStubRef->getDetId() gives the stackDetId, not rawId
      DetId detIdStub = theTrackerGeometry->idToDet( (tempStubRef->getClusterRef(0))->getDetId() )->geographicalId();
      
      // CHECK IF THIS STILL WORKS !!
      bool genuineStub    = MCTruthTTStubHandle->isGenuine( tempStubRef );
      bool combinStub     = MCTruthTTStubHandle->isCombinatoric( tempStubRef );
      
      /// Define position stub by position inner cluster
      MeasurementPoint mp = (tempStubRef->getClusterRef(0))->findAverageLocalCoordinates();
      const GeomDet* theGeomDet = theTrackerGeometry->idToDet(detIdStub);
      Global3DPoint posStub = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition(mp) );
      
      if ( detIdStub.subdetId() == static_cast<int>(StripSubdetector::TOB) )  // Phase 2 Outer Tracker Barrel
      {
        if ( genuineStub )
        {
          Stub_Gen_Barrel->Fill( tTopo->layer(detIdStub) );
        }
        else if ( combinStub )
        {
          Stub_Comb_Barrel->Fill( tTopo->layer(detIdStub) );
        }
        else
        {
          Stub_Unkn_Barrel->Fill( tTopo->layer(detIdStub) );
        } 
      } // end if isBarrel
      else if ( detIdStub.subdetId() == static_cast<int>(StripSubdetector::TID) )  // Phase 2 Outer Tracker Endcap
      {
        int side = tTopo->side(detIdStub);
        int disc = tTopo->layer(detIdStub);  // returns wheel
        int ring = tTopo->tidRing(detIdStub);
        if ( genuineStub )
        {
          Stub_Gen_Endcap_Disc->Fill( disc );
          Stub_Gen_Endcap_Ring->Fill( ring );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Stub_Gen_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Stub_Gen_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
        else if ( combinStub )
        {
          Stub_Comb_Endcap_Disc->Fill( disc );
          Stub_Comb_Endcap_Ring->Fill( ring );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Stub_Comb_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Stub_Comb_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
        else
        {
          Stub_Unkn_Endcap_Disc->Fill( disc );
          Stub_Unkn_Endcap_Ring->Fill( ring );
          if ( verbosePlots_ )
          {
            if ( side == 1 ) Stub_Unkn_Endcap_Ring_Bw[disc-1]->Fill( ring );
            else if ( side == 2 ) Stub_Unkn_Endcap_Ring_Fw[disc-1]->Fill( ring );
          }  /// end verbosePlots
        }
      }	// end if isEndcap
      
      /// Eta distribution in function of genuine/combinatorial/unknown stub
      if ( genuineStub ) Stub_Gen_Eta->Fill( posStub.eta() );
      else if ( combinStub ) Stub_Comb_Eta->Fill( posStub.eta() );
      else Stub_Unkn_Eta->Fill( posStub.eta() );
      
    }	// end loop contentIter
  }	// end loop inputIter
}


// ------------ method called when starting to processes a run  ------------

void 
OuterTrackerStub::beginRun(edm::Run const&, edm::EventSetup const&)
{
  std::string HistoName;
  
  dqmStore_->setCurrentFolder(topFolderName_+"/TTStubs/");
  
  edm::ParameterSet psTTStubEta =  conf_.getParameter<edm::ParameterSet>("TH1TTStub_Eta");
  HistoName = "Stub_Gen_Eta";
  Stub_Gen_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTStubEta.getParameter<int32_t>("Nbinsx"),
      psTTStubEta.getParameter<double>("xmin"),
      psTTStubEta.getParameter<double>("xmax"));
  Stub_Gen_Eta->setAxisTitle("#eta", 1);
  Stub_Gen_Eta->setAxisTitle("# Genuine L1 Stubs", 2);
  
  HistoName = "Stub_Unkn_Eta";
  Stub_Unkn_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTStubEta.getParameter<int32_t>("Nbinsx"),
      psTTStubEta.getParameter<double>("xmin"),
      psTTStubEta.getParameter<double>("xmax"));
  Stub_Unkn_Eta->setAxisTitle("#eta", 1);
  Stub_Unkn_Eta->setAxisTitle("# Unknown L1 Stubs", 2);
  
  HistoName = "Stub_Comb_Eta";
  Stub_Comb_Eta = dqmStore_->book1D(HistoName, HistoName,
      psTTStubEta.getParameter<int32_t>("Nbinsx"),
      psTTStubEta.getParameter<double>("xmin"),
      psTTStubEta.getParameter<double>("xmax"));
  Stub_Comb_Eta->setAxisTitle("#eta", 1);
  Stub_Comb_Eta->setAxisTitle("# Combinatorial L1 Stubs", 2);
  
  
  /// TTStub stacks
  edm::ParameterSet psTTStubLayer =  conf_.getParameter<edm::ParameterSet>("TH1TTStub_Layer");
  HistoName = "NStubs_Gen_Barrel";
  Stub_Gen_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTStubLayer.getParameter<int32_t>("Nbinsx"),
      psTTStubLayer.getParameter<double>("xmin"),
      psTTStubLayer.getParameter<double>("xmax"));
  Stub_Gen_Barrel->setAxisTitle("Barrel Layer", 1);
  Stub_Gen_Barrel->setAxisTitle("# Genuine L1 Stubs", 2);
  
  HistoName = "NStubs_Unkn_Barrel";
  Stub_Unkn_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTStubLayer.getParameter<int32_t>("Nbinsx"),
      psTTStubLayer.getParameter<double>("xmin"),
      psTTStubLayer.getParameter<double>("xmax"));
  Stub_Unkn_Barrel->setAxisTitle("Barrel Layer", 1);
  Stub_Unkn_Barrel->setAxisTitle("# Unknown L1 Stubs", 2);
  
  HistoName = "NStubs_Comb_Barrel";
  Stub_Comb_Barrel = dqmStore_->book1D(HistoName, HistoName,
      psTTStubLayer.getParameter<int32_t>("Nbinsx"),
      psTTStubLayer.getParameter<double>("xmin"),
      psTTStubLayer.getParameter<double>("xmax"));
  Stub_Comb_Barrel->setAxisTitle("Barrel Layer", 1);
  Stub_Comb_Barrel->setAxisTitle("# Combinatorial L1 Stubs", 2);
  
  edm::ParameterSet psTTStubDisk =  conf_.getParameter<edm::ParameterSet>("TH1TTStub_Disk");
  HistoName = "NStubs_Gen_Endcap_Disc";
  Stub_Gen_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTStubDisk.getParameter<int32_t>("Nbinsx"),
      psTTStubDisk.getParameter<double>("xmin"),
      psTTStubDisk.getParameter<double>("xmax"));
  Stub_Gen_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Stub_Gen_Endcap_Disc->setAxisTitle("# Genuine L1 Stubs", 2);
  
  HistoName = "NStubs_Unkn_Endcap_Disc";
  Stub_Unkn_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTStubDisk.getParameter<int32_t>("Nbinsx"),
      psTTStubDisk.getParameter<double>("xmin"),
      psTTStubDisk.getParameter<double>("xmax"));
  Stub_Unkn_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Stub_Unkn_Endcap_Disc->setAxisTitle("# Unknown L1 Stubs", 2);
  
  HistoName = "NStubs_Comb_Endcap_Disc";
  Stub_Comb_Endcap_Disc = dqmStore_->book1D(HistoName, HistoName,
      psTTStubDisk.getParameter<int32_t>("Nbinsx"),
      psTTStubDisk.getParameter<double>("xmin"),
      psTTStubDisk.getParameter<double>("xmax"));
  Stub_Comb_Endcap_Disc->setAxisTitle("Endcap Disc", 1);
  Stub_Comb_Endcap_Disc->setAxisTitle("# Combinatorial L1 Stubs", 2);
  
  edm::ParameterSet psTTStubRing =  conf_.getParameter<edm::ParameterSet>("TH1TTStub_Ring");
  HistoName = "NStubs_Gen_Endcap_Ring";
  Stub_Gen_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTStubRing.getParameter<int32_t>("Nbinsx"),
      psTTStubRing.getParameter<double>("xmin"),
      psTTStubRing.getParameter<double>("xmax"));
  Stub_Gen_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Stub_Gen_Endcap_Ring->setAxisTitle("# Genuine L1 Stubs", 2);
  
  HistoName = "NStubs_Unkn_Endcap_Ring";
  Stub_Unkn_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTStubRing.getParameter<int32_t>("Nbinsx"),
      psTTStubRing.getParameter<double>("xmin"),
      psTTStubRing.getParameter<double>("xmax"));
  Stub_Unkn_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Stub_Unkn_Endcap_Ring->setAxisTitle("# Unknown L1 Stubs", 2);
  
  HistoName = "NStubs_Comb_Endcap_Ring";
  Stub_Comb_Endcap_Ring = dqmStore_->book1D(HistoName, HistoName,
      psTTStubRing.getParameter<int32_t>("Nbinsx"),
      psTTStubRing.getParameter<double>("xmin"),
      psTTStubRing.getParameter<double>("xmax"));
  Stub_Comb_Endcap_Ring->setAxisTitle("Endcap Ring", 1);
  Stub_Comb_Endcap_Ring->setAxisTitle("# Combinatorial L1 Stubs", 2);
  
  
  /// Plots for debugging
  if ( verbosePlots_ )
  {
    dqmStore_->setCurrentFolder(topFolderName_+"/TTStubs/NStubsPerRing");

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Gen_Disk+"+std::to_string(i+1);
      Stub_Gen_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Gen_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Gen_Endcap_Ring_Fw[i]->setAxisTitle("# Genuine L1 Stubs", 2);
    }  

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Gen_Disk-"+std::to_string(i+1);
      Stub_Gen_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Gen_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Gen_Endcap_Ring_Bw[i]->setAxisTitle("# Genuine L1 Stubs", 2);
    }

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Unkn_Disk+"+std::to_string(i+1);
      Stub_Unkn_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("# Unknown L1 Stubs", 2);
    }

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Unkn_Disk-"+std::to_string(i+1);
      Stub_Unkn_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("# Unknown L1 Stubs", 2);
    }

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Comb_Disk+"+std::to_string(i+1);
      Stub_Comb_Endcap_Ring_Fw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Comb_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Comb_Endcap_Ring_Fw[i]->setAxisTitle("# Combinatorial L1 Stubs", 2);
    }

    for (int i = 0; i < 5; i++)
    {
      HistoName = "NStubs_Comb_Disk-"+std::to_string(i+1);
      Stub_Comb_Endcap_Ring_Bw[i] = dqmStore_->book1D(HistoName, HistoName,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Comb_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Comb_Endcap_Ring_Bw[i]->setAxisTitle("# Combinatorial L1 Stubs", 2);
    }
    
  }  /// end verbosePlots
}


// ------------ method called once each job just after ending the event loop  ------------
void 
OuterTrackerStub::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(OuterTrackerStub);
