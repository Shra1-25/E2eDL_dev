#ifndef DetFrameProducer_h
#define DetFrameProducer_h

#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"

#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TFrame.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TVector3.h"
#include "TVector2.h"
#include <iostream>

#include "E2eDL/DataFormats/interface/FrameCollections.h"
//#include "E2eDL/FrameProducers/interface/constants.h"
#include "E2eDL/FrameProducers/interface/fillECALstitched.h"
#include "E2eDL/FrameProducers/interface/fillHBHE.h"
#include "E2eDL/FrameProducers/interface/fillTracksAtECALadj.h"

using namespace std;
/*using pat::PhotonCollection;
using pat::PhotonRef;*/
using reco::PhotonCollection;
using reco::PhotonRef;

class DetFrameProducer : public edm::stream::EDProducer<> {
   public:
      
      explicit DetFrameProducer(const edm::ParameterSet&);
      ~DetFrameProducer();
      
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      // Tokens
      edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_; 
      edm::EDGetTokenT<PhotonCollection> photonCollectionT_;
      edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
      edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
      edm::EDGetTokenT<reco::TrackCollection> trackCollectionT_;
      edm::EDGetTokenT<reco::VertexCollection> vertexCollectionT_;
      edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
      edm::EDGetTokenT<reco::JetTagCollection> jetTagCollectionT_;
      edm::EDGetTokenT<std::vector<reco::CandIPTagInfo> >    ipTagInfoCollectionT_;
      edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
      edm::EDGetTokenT<TrackingRecHitCollection> TRKRecHitCollectionT_;
      edm::EDGetTokenT<edm::View<reco::Jet> > recoJetsT_;
      edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
   
      // Detector image switches
      bool doECALstitched;
      bool doTracksAtECALstitchedPt;
      bool doTracksAtECALadjPt;
      bool doHBHEenergy;
      
      static const int nPhotons = 2;
      
      TProfile2D *hEB_energy;
      TProfile2D *hEB_time;
      //TProfile2D *hEB_frame;
      e2e::Frame1D vEB_energy_;
      e2e::Frame1D vEB_time_;
      e2e::Frame1D vHBHE_energy_EB_;
      e2e::Frame1D vHBHE_energy_;
      e2e::Frame2D vHBHE_energy_reshaped;
      e2e::Frame1D vECAL_energy_;
      e2e::Frame2D vECAL_energy_reshaped; 
      e2e::Frame1D vECAL_tracksPt_;
      e2e::Frame2D vECAL_tracksPt_reshaped;
      e2e::Frame1D vECALadj_tracksPt_[Nadjproj];
      e2e::Frame2D vECALadj_tracksPt_reshaped;
      e2e::Frame1D vECALadj_tracks_[Nadjproj];
      e2e::Frame1D vECALadj_tracksPt_max_[Nadjproj];
      
      unsigned int nPho;
      
      /*typedef reco::VertexCollection  PVCollection;
      edm::EDGetTokenT<PVCollection> pvCollectionT_;*/
      typedef std::vector<reco::PFCandidate>  PFCollection;
      edm::EDGetTokenT<PFCollection> pfCollectionT_;
      
      void fillEB             ( const edm::Event&, const edm::EventSetup& );
      void fillTracksAtECALstitched (const edm::Event&, const edm::EventSetup& );
      void fillTracksAtECALadjustable   ( const edm::Event&, const edm::EventSetup&, unsigned int proj );
      
      std::vector<int> findSubcrystal(const CaloGeometry* caloGeom, const float& eta, const float& phi, const int& granularityMultiEta, const int& granularityMultiPhi);
      void fillByBinNumber(TH2F * histo, const std::vector<int>& phi_eta, const float& value);   
      
      std::vector<float>& read_vEB_energy     (int);
      
      int iphi_Emax, ieta_Emax;
      /*unsigned int granularityMultiPhi[Nadjproj];
      unsigned int granularityMultiEta[Nadjproj];
      
      int totalEtaBins[Nadjproj];// = totalMultiEta*(eta_nbins_HBHE);
      int totalPhiBins[Nadjproj];// = granularityMultiPhi * granularityMultiECAL*HBHE_IPHI_NUM;
      std::vector<double> adjEtaBins[Nadjproj];*/
      
      e2e::Frame1D vIphi_Emax_;
      e2e::Frame1D vIeta_Emax_;
      
      std::vector<int> vPreselPhoIdxs_;
      int nTotal=0;
      int nPassed=0;
};

#endif
//DEFINE_FWK_MODULE(DetFrameProducer);
