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
#include "E2eDL/FrameProducers/interface/constants.h"
#include "E2eDL/FrameProducers/interface/fillECALstitched.h"
#include "E2eDL/FrameProducers/interface/fillHBHE.h"
#include "E2eDL/FrameProducers/interface/fillTracksAtECALadj.h"

using namespace std;
/*using pat::PhotonCollection;
using pat::PhotonRef;*/
using reco::PhotonCollection;
using reco::PhotonRef;

static const unsigned int Nproj = 5;
//static const unsigned int Nadjproj = 2;
static const unsigned int Nhitproj = 2;
static const unsigned int nDetFrameH = 280;
static const unsigned int nDetFrameW = 360;

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
      e2e::Frame2D vECALad_tracksPt_reshaped;
      e2e::Frame1D vECALadj_tracks_[Nadjproj];
      e2e::Frame1D vECALadj_tracksPt_max_[Nadjproj];
      
      unsigned int nPho;
      
      /*typedef reco::VertexCollection  PVCollection;
      edm::EDGetTokenT<PVCollection> pvCollectionT_;*/
      typedef std::vector<reco::PFCandidate>  PFCollection;
      edm::EDGetTokenT<PFCollection> pfCollectionT_;
      
      void fillEB             ( const edm::Event&, const edm::EventSetup& );
      void fillTracksAtECALstitched (const edm::Event&, const edm::EventSetup& );
      //void fillTracksAtECALadjustable   ( const edm::Event&, const edm::EventSetup&, unsigned int proj );
      
      //std::vector<int> findSubcrystal(const CaloGeometry* caloGeom, const float& eta, const float& phi, const int& granularityMultiEta, const int& granularityMultiPhi);
      //void fillByBinNumber(TH2F * histo, const std::vector<int>& phi_eta, const float& value);   
      
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

static const bool debug = false;


static const int nTOB = 6;
static const int nTEC = 9;
static const int nTIB = 4;
static const int nTID = 3;
static const int nBPIX = 4;
static const int nFPIX = 3;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100
static const int HBHE_IETA_MAX_FINE = 20;
static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;
static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1;
static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;


static const std::string projections[Nproj] = {"", "_atECAL", "_atHCAL","_atECALfixIP","_atECALfixIPfromPV"}; //57425
static const std::string hit_projections[Nhitproj] = {"", "_atPV"};
static const int eta_nbins_HBHE = 2*(HBHE_IETA_MAX_HE-1);
static const int granularityMultiECAL=5;

// HBHE eta bin edges
static const double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
                  {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                   -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57
//
// static data member definitions
//

#endif
//DEFINE_FWK_MODULE(DetFrameProducer);
