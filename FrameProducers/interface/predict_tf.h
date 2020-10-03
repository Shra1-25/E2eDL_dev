#ifndef RECO_E2E_predict_tf_h
#define RECO_E2E_predict_tf_h

#include <vector>
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "E2eDL/DataFormats/interface/FrameCollections.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "tensorflow/core/graph/default_device.h"
#include <iostream>
#include <fstream>
using namespace std;

namespace e2e {
  e2e::Frame2D predict_tf(e2e::Frame4D&, string, string, string);
}
#endif
