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

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
#include "Validation/Phase2OuterTracker/interface/OuterTrackerStub.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"	//Needed for TTStubAssociationMap.h !
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "TMath.h"
#include <iostream>

//
// constructors and destructor
//
OuterTrackerStub::OuterTrackerStub(const edm::ParameterSet& iConfig)
: dqmStore_(edm::Service<DQMStore>().operator->()), conf_(iConfig)

{
  topFolderName_ = conf_.getParameter<std::string>("TopFolderName");
  tagTTStubs_ = conf_.getParameter< edm::InputTag >("TTStubs");
  tagTTStubMCTruth_ = conf_.getParameter< edm::InputTag >("TTStubMCTruth");
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
  /// Geometry handles etc
  edm::ESHandle< TrackerGeometry >                GeometryHandle;
  edm::ESHandle< StackedTrackerGeometry >         StackedGeometryHandle;
  const StackedTrackerGeometry*                   theStackedGeometry;
  
  /// Geometry setup
  /// Set pointers to Geometry
  iSetup.get< TrackerDigiGeometryRecord >().get(GeometryHandle);
  /// Set pointers to Stacked Modules
  iSetup.get< StackedTrackerGeometryRecord >().get(StackedGeometryHandle);
  theStackedGeometry = StackedGeometryHandle.product(); /// Note this is different from the "global" geometry
  
  /// Track Trigger
  edm::Handle< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > > > PixelDigiTTStubHandle;
  iEvent.getByLabel( tagTTStubs_, PixelDigiTTStubHandle );
  /// Track Trigger MC Truth
  edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > > MCTruthTTStubHandle;
  iEvent.getByLabel( tagTTStubMCTruth_, MCTruthTTStubHandle );
  
  /// Loop over the input Stubs
  typename edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >::const_iterator otherInputIter;
  typename edmNew::DetSet< TTStub< Ref_PixelDigi_ > >::const_iterator otherContentIter;
  for ( otherInputIter = PixelDigiTTStubHandle->begin();
       otherInputIter != PixelDigiTTStubHandle->end();
       ++otherInputIter )
  {
    for ( otherContentIter = otherInputIter->begin();
         otherContentIter != otherInputIter->end();
         ++otherContentIter )
    {
      /// Make the reference to be put in the map
      edm::Ref< edmNew::DetSetVector< TTStub< Ref_PixelDigi_ > >, TTStub< Ref_PixelDigi_ > > tempStubRef = edmNew::makeRefTo( PixelDigiTTStubHandle, otherContentIter );
      
      StackedTrackerDetId detIdStub( tempStubRef->getDetId() );
      
      bool genuineStub    = MCTruthTTStubHandle->isGenuine( tempStubRef );
      bool combinStub     = MCTruthTTStubHandle->isCombinatoric( tempStubRef );
      
      GlobalPoint posStub = theStackedGeometry->findGlobalPosition( &(*tempStubRef) );
      
      if ( detIdStub.isBarrel() )
      {
        if ( genuineStub )
        {
          Stub_Gen_Barrel->Fill( detIdStub.iLayer() );
        }
        else if ( combinStub )
        {
          Stub_Comb_Barrel->Fill( detIdStub.iLayer() );
        }
        else
        {
          Stub_Unkn_Barrel->Fill( detIdStub.iLayer() );
        } 
      } // end if isBarrel()
      else if ( detIdStub.isEndcap() )
      {
        if ( genuineStub )
        {
          Stub_Gen_Endcap_Disc->Fill( detIdStub.iDisk() );
          Stub_Gen_Endcap_Ring->Fill( detIdStub.iRing() );
          if ( verbosePlots_ )
          {
            if ( detIdStub.iSide() == 1) Stub_Gen_Endcap_Ring_Bw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
            else if ( detIdStub.iSide() == 2) Stub_Gen_Endcap_Ring_Fw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
          } /// End verbosePlots
        }
        else if ( combinStub )
        {
          Stub_Comb_Endcap_Disc->Fill( detIdStub.iDisk() );
          Stub_Comb_Endcap_Ring->Fill( detIdStub.iRing() );
          if ( verbosePlots_ )
          {
            if ( detIdStub.iSide() == 1) Stub_Comb_Endcap_Ring_Bw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
            else if ( detIdStub.iSide() == 2) Stub_Comb_Endcap_Ring_Fw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
          } /// End verbosePlots
        }
        else
        {
          Stub_Unkn_Endcap_Disc->Fill( detIdStub.iDisk() );
          Stub_Unkn_Endcap_Ring->Fill( detIdStub.iRing() );
          if ( verbosePlots_ )
          {
            if ( detIdStub.iSide() == 1) Stub_Unkn_Endcap_Ring_Bw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
            else if ( detIdStub.iSide() == 2) Stub_Unkn_Endcap_Ring_Fw[detIdStub.iDisk()-1]->Fill( detIdStub.iRing() );
          } /// End verbosePlots
        }
      }	// end if isEndcap()
      
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
  SiStripFolderOrganizer folder_organizer;
  folder_organizer.setSiStripFolderName(topFolderName_);
  folder_organizer.setSiStripFolder();
  std::string HistoName;
  
  dqmStore_->setCurrentFolder(topFolderName_+"/Stubs/");
  
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
    dqmStore_->setCurrentFolder(topFolderName_+"/Stubs/NStubsPerRing");

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Gen_Disk+%d", i+1);
      Stub_Gen_Endcap_Ring_Fw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Gen_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Gen_Endcap_Ring_Fw[i]->setAxisTitle("# Genuine L1 Stubs", 2);
    }  

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Gen_Disk-%d", i+1);
      Stub_Gen_Endcap_Ring_Bw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Gen_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Gen_Endcap_Ring_Bw[i]->setAxisTitle("# Genuine L1 Stubs", 2);
    }

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Unkn_Disk+%d", i+1);
      Stub_Unkn_Endcap_Ring_Fw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Unkn_Endcap_Ring_Fw[i]->setAxisTitle("# Unknown L1 Stubs", 2);
    }

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Unkn_Disk-%d", i+1);
      Stub_Unkn_Endcap_Ring_Bw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Unkn_Endcap_Ring_Bw[i]->setAxisTitle("# Unknown L1 Stubs", 2);
    }

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Comb_Disk+%d", i+1);
      Stub_Comb_Endcap_Ring_Fw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Comb_Endcap_Ring_Fw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Comb_Endcap_Ring_Fw[i]->setAxisTitle("# Combinatorial L1 Stubs", 2);
    }

    for(int i=0;i<5;i++){
      Char_t histo[200];
      sprintf(histo, "NStubs_Comb_Disk-%d", i+1);
      Stub_Comb_Endcap_Ring_Bw[i] = dqmStore_->book1D(histo, histo,
          psTTStubRing.getParameter<int32_t>("Nbinsx"),
          psTTStubRing.getParameter<double>("xmin"),
          psTTStubRing.getParameter<double>("xmax"));
      Stub_Comb_Endcap_Ring_Bw[i]->setAxisTitle("Endcap Ring", 1);
      Stub_Comb_Endcap_Ring_Bw[i]->setAxisTitle("# Combinatorial L1 Stubs", 2);
    }
    
  } /// End verbosePlots
}


// ------------ method called once each job just after ending the event loop  ------------
void 
OuterTrackerStub::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(OuterTrackerStub);
