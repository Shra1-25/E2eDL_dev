#ifndef RecoE2E_fillECALstitched_h
#define RecoE2E_fillECALstitched_h

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
//#include "E2eDL/FrameProducers/interface/DetFrameProducer.h"

namespace e2e {
  e2e::Frame1D fillECALstitched   ( edm::Handle<EcalRecHitCollection> EBRecHitsH_, edm::Handle<EcalRecHitCollection> EERecHitsH_, edm::ESHandle<CaloGeometry> caloGeomH_ );
}
static const int nEE = 2;
static const float zs = 0.;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const int ECAL_IETA_MAX_EXT = 140;
static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
static const double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                  {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                   -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                   -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                   -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                   -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                   -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                   -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                   -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
// EE+(phi,eta) projection eta edges
static const double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                   {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                    1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                    1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                    1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                    1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                    2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                    2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                    2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56

#endif
