#ifndef Phase2OuterTracker_OuterTrackerMCTruth_h
#define Phase2OuterTracker_OuterTrackerMCTruth_h

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
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
//#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"  // does not exist yet in 90X
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
//#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"  // does not exist yet in 90X
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

class PixelGeomDetUnit;
class PixelTopology;
using Phase2TrackerGeomDetUnit = PixelGeomDetUnit;
using Phase2TrackerTopology = PixelTopology;

class DQMStore;

class OuterTrackerMCTruth : public edm::EDAnalyzer {

public:
  explicit OuterTrackerMCTruth(const edm::ParameterSet&);
  ~OuterTrackerMCTruth();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&);

  template< typename T > LocalPoint findClusterLocalPosition(const TTCluster<T> *cluster, const TrackerGeometry *tkGeom);
  template< typename T > GlobalPoint findStubGlobalPosition(const TTStub<T> *stub, const TrackerGeometry *tkGeom, unsigned int stackMember);
  template< typename T > GlobalVector findStubGlobalDirection(const TTStub<T> *stub, const TrackerGeometry *tkGeom);
  template< typename T > double getRoughStubPt(const TTStub<T> *stub, const TrackerGeometry *tkGeom, double magnFieldStrength);
 
  
  // TrackingParticle and TrackingVertex
  MonitorElement* SimVtx_XY = 0;
  MonitorElement* SimVtx_RZ = 0;
  
  MonitorElement* TPart_Pt = 0;
  MonitorElement* TPart_Eta_Pt10 = 0;
  MonitorElement* TPart_Phi_Pt10 = 0;
  
  MonitorElement* TPart_Cluster_Pt = 0;
  MonitorElement* TPart_Cluster_Phi_Pt10 = 0;
  MonitorElement* TPart_Cluster_Eta_Pt10 = 0;
  
  MonitorElement* TPart_Stub_Pt = 0;
  MonitorElement* TPart_Stub_Phi_Pt10 = 0;
  MonitorElement* TPart_Stub_Eta_Pt10 = 0;
  
  MonitorElement* TPart_Track_LQ_Pt = 0;
  MonitorElement* TPart_Track_LQ_Phi_Pt10 = 0;
  MonitorElement* TPart_Track_LQ_Eta_Pt10 = 0;
  MonitorElement* TPart_Track_HQ_Pt = 0;
  MonitorElement* TPart_Track_HQ_Phi_Pt10 = 0;
  MonitorElement* TPart_Track_HQ_Eta_Pt10 = 0;
  
  // Stub in PS/2S module vs. TPart Eta
  MonitorElement* TPart_Stub_Eta_Pt10_Normalization = 0;
  MonitorElement* TPart_Stub_Eta_Pt10_NumPS = 0;
  MonitorElement* TPart_Stub_Eta_Pt10_Num2S = 0;
  
  // CW vs. TPart Eta
  MonitorElement* TPart_Eta_INormalization = 0;
  MonitorElement* TPart_Eta_ICW_1 = 0;
  MonitorElement* TPart_Eta_ICW_2 = 0;
  MonitorElement* TPart_Eta_ICW_3 = 0;
  MonitorElement* TPart_Eta_ONormalization = 0;
  MonitorElement* TPart_Eta_OCW_1 = 0;
  MonitorElement* TPart_Eta_OCW_2 = 0;
  MonitorElement* TPart_Eta_OCW_3 = 0;
  
  // PID
  MonitorElement* Cluster_PID = 0;
  MonitorElement* Stub_PID = 0;
  
  // Track Chi2(Red)
  MonitorElement* Track_LQ_Chi2_TPart_Eta = 0;
  MonitorElement* Track_LQ_Chi2Red_TPart_Eta = 0;
  MonitorElement* Track_HQ_Chi2_TPart_Eta = 0;
  MonitorElement* Track_HQ_Chi2Red_TPart_Eta = 0;
  
  
  // Stubs vs. TrackingParticles
  MonitorElement* Stub_InvPt_TPart_InvPt_AllLayers = 0;
  MonitorElement* Stub_Pt_TPart_Pt_AllLayers = 0;
  MonitorElement* Stub_Eta_TPart_Eta_AllLayers = 0;
  MonitorElement* Stub_Phi_TPart_Phi_AllLayers = 0;
  
  MonitorElement* Stub_InvPtRes_TPart_Eta_AllLayers = 0;
  MonitorElement* Stub_PtRes_TPart_Eta_AllLayers = 0;
  MonitorElement* Stub_EtaRes_TPart_Eta_AllLayers = 0;
  MonitorElement* Stub_PhiRes_TPart_Eta_AllLayers = 0;
  MonitorElement* Stub_PhiRes_TPart_Phi_AllLayers = 0;
  
  MonitorElement* Stub_W_TPart_Pt_AllLayers = 0;
  MonitorElement* Stub_W_TPart_InvPt_AllLayers = 0;
  
  MonitorElement* Stub_InvPt_TPart_InvPt_AllDiscs = 0;
  MonitorElement* Stub_Pt_TPart_Pt_AllDiscs = 0;
  MonitorElement* Stub_Eta_TPart_Eta_AllDiscs = 0;
  MonitorElement* Stub_Phi_TPart_Phi_AllDiscs = 0;
  
  MonitorElement* Stub_InvPtRes_TPart_Eta_AllDiscs = 0;
  MonitorElement* Stub_PtRes_TPart_Eta_AllDiscs = 0;
  MonitorElement* Stub_EtaRes_TPart_Eta_AllDiscs = 0;
  MonitorElement* Stub_PhiRes_TPart_Eta_AllDiscs = 0;
  MonitorElement* Stub_PhiRes_TPart_Phi_AllDiscs = 0;
  
  MonitorElement* Stub_W_TPart_Pt_AllDiscs = 0;
  MonitorElement* Stub_W_TPart_InvPt_AllDiscs = 0;
  
  
  // Tracks vs. TrackingParticles
  MonitorElement* Track_LQ_Pt_TPart_Pt = 0;
  MonitorElement* Track_LQ_PtRes_TPart_Eta = 0;
  MonitorElement* Track_LQ_InvPt_TPart_InvPt = 0;
  MonitorElement* Track_LQ_InvPtRes_TPart_Eta = 0;
  MonitorElement* Track_LQ_Phi_TPart_Phi = 0;
  MonitorElement* Track_LQ_PhiRes_TPart_Eta = 0;
  MonitorElement* Track_LQ_Eta_TPart_Eta = 0;
  MonitorElement* Track_LQ_EtaRes_TPart_Eta = 0;
  MonitorElement* Track_LQ_VtxZ0_TPart_VtxZ0 = 0;
  MonitorElement* Track_LQ_VtxZ0Res_TPart_Eta = 0;
  
  MonitorElement* Track_HQ_Pt_TPart_Pt = 0;
  MonitorElement* Track_HQ_PtRes_TPart_Eta = 0;
  MonitorElement* Track_HQ_InvPt_TPart_InvPt = 0;
  MonitorElement* Track_HQ_InvPtRes_TPart_Eta = 0;
  MonitorElement* Track_HQ_Phi_TPart_Phi = 0;
  MonitorElement* Track_HQ_PhiRes_TPart_Eta = 0;
  MonitorElement* Track_HQ_Eta_TPart_Eta = 0;
  MonitorElement* Track_HQ_EtaRes_TPart_Eta = 0;
  MonitorElement* Track_HQ_VtxZ0_TPart_VtxZ0 = 0;
  MonitorElement* Track_HQ_VtxZ0Res_TPart_Eta = 0;
  
  
 private:
  DQMStore* dqmStore_;
  edm::ParameterSet conf_;
  edm::EDGetTokenT< TrackingParticleCollection > tagTParticleToken_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTCluster< Ref_Phase2TrackerDigi_ > > >  tagTTClustersToken_;
  edm::EDGetTokenT< TTClusterAssociationMap< Ref_Phase2TrackerDigi_ > >  tagTTClusterMCTruthToken_;
  edm::EDGetTokenT<edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >  tagTTStubsToken_;
  edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >  tagTTStubMCTruthToken_;
  //edm::EDGetTokenT<edmNew::DetSetVector< TTTrack< Ref_Phase2TrackerDigi_ > > >  tagTTTracksToken_;
  //edm::EDGetTokenT< TTTrackAssociationMap< Ref_Phase2TrackerDigi_ > >  tagTTTrackMCTruthToken_;

  std::string topFolderName_;
  unsigned int HQDelim_;
  bool verbosePlots_;
  
};

template< typename T >
LocalPoint OuterTrackerMCTruth::findClusterLocalPosition(const TTCluster<T> *cluster, const TrackerGeometry *tkGeom)
{
  const GeomDet* theGeomDet = tkGeom->idToDet( cluster->getDetId() );
  return theGeomDet->topology().localPosition( cluster->findAverageLocalCoordinates() );
}

template< typename T >
GlobalPoint OuterTrackerMCTruth::findStubGlobalPosition(const TTStub<T> *stub, const TrackerGeometry *tkGeom, unsigned int stackMember)
{
  const GeomDet* theGeomDet = tkGeom->idToDet( (stub->getClusterRef(stackMember))->getDetId() );
  const GlobalPoint thePosition = theGeomDet->surface().toGlobal( theGeomDet->topology().localPosition( (stub->getClusterRef(stackMember))->findAverageLocalCoordinates() ) );

  return thePosition;
}

template< typename T >
GlobalVector OuterTrackerMCTruth::findStubGlobalDirection(const TTStub<T> *stub, const TrackerGeometry *tkGeom)
{
  const GlobalPoint innerHitPosition = findStubGlobalPosition(stub, tkGeom, 0);
  const GlobalPoint outerHitPosition = findStubGlobalPosition(stub, tkGeom, 1);
  const GlobalVector tempStubDirection( outerHitPosition.x()-innerHitPosition.x(), outerHitPosition.y()-innerHitPosition.y(), outerHitPosition.z()-innerHitPosition.z() );

  return tempStubDirection;
}

template< typename T > double OuterTrackerMCTruth::getRoughStubPt(const TTStub<T> *stub, const TrackerGeometry *tkGeom, double magnFieldStrength )
{
  /// Get the magnetic field
  double mPtFactor = (floor(magnFieldStrength*10.0 + 0.5))/10.0*0.0015;;  // CHECK !!

  /// Get average position of Clusters composing the Stub
  const GlobalPoint innerHitPosition = findStubGlobalPosition(stub, tkGeom, 0);
  const GlobalPoint outerHitPosition = findStubGlobalPosition(stub, tkGeom, 1);

  /// Get useful quantities
  double outerPointRadius = outerHitPosition.perp();
  double innerPointRadius = innerHitPosition.perp();
  double deltaRadius = outerPointRadius - innerPointRadius;

  DetId tempDetId = tkGeom->idToDet( (stub->getClusterRef(0))->getDetId() )->geographicalId();
  if ( tempDetId.subdetId() == static_cast<int>(StripSubdetector::TOB) )  // Phase 2 Outer Tracker Barrel
  {
    /// Calculate angular displacement from hit phi locations
    /// and renormalize it, if needed
    double deltaPhi = outerHitPosition.phi() - innerHitPosition.phi();
    if (deltaPhi < 0)
      deltaPhi = -deltaPhi;
    else if (deltaPhi > M_PI)
      deltaPhi = 2*M_PI - deltaPhi;

    /// Return the rough Pt
    return ( deltaRadius * mPtFactor / deltaPhi );
  }
  else if ( tempDetId.subdetId() == static_cast<int>(StripSubdetector::TID) )  // Phase 2 Outer Tracker Endcap
  {
    /// Test approximated formula for Endcap stubs
    /// Check always to be consistent with HitMatchingAlgorithm_window2012.h
    double roughPt = innerPointRadius * innerPointRadius * mPtFactor / fabs( findClusterLocalPosition(&(*stub->getClusterRef(0)), tkGeom).x() );
    roughPt +=       outerPointRadius * outerPointRadius * mPtFactor / fabs( findClusterLocalPosition(&(*stub->getClusterRef(1)), tkGeom).x() );
    roughPt = roughPt / 2.;

    return roughPt;
  }

  /// Should never come here (! barrel && ! endcap)
  return -9999.;
}
#endif
