#ifndef RecoE2E_fillECALstitched_h
#define RecoE2E_fillECALstitched_h
#include "E2eDL/FrameProducers/interface/DetFrameProducer.h"

namespace e2e {
  e2e::Frame1D fillECALstitched   ( edm::Handle<EcalRecHitCollection> EBRecHitsH_, edm::Handle<EcalRecHitCollection> EERecHitsH_, edm::ESHandle<CaloGeometry> caloGeomH_ );
}
#endif
