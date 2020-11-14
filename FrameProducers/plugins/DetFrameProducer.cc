//DetFrameProducer
#include "E2eDL/FrameProducers/interface/DetFrameProducer.h"


DetFrameProducer::DetFrameProducer(const edm::ParameterSet& iConfig)
{
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  photonCollectionT_ = consumes<PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonCollection"));
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
  EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  trackCollectionT_       = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  vertexCollectionT_       = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackRecHitCollection"));
  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  pfCollectionT_          = consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCollection"));
  pvCollectionT_ = consumes<PVCollection>(iConfig.getParameter<edm::InputTag>("pvCollection"));
  jetTagCollectionT_      = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("jetTagCollection"));
  ipTagInfoCollectionT_   = consumes<std::vector<reco::CandIPTagInfo> > (iConfig.getParameter<edm::InputTag>("ipTagInfoCollection"));
  
  ///////////adjustable granularity stuff

  granularityMultiPhi[0]  = iConfig.getParameter<int>("granularityMultiPhi");
  granularityMultiEta[0]  = iConfig.getParameter<int>("granularityMultiEta");

  granularityMultiPhi[1] = 3;
  granularityMultiEta[1] = 3;
  
  for (unsigned int proj=0; proj<Nadjproj; proj++)
  {
  
    int totalMultiEta = granularityMultiEta[proj] * granularityMultiECAL;

    for (int i=0; i<eta_nbins_HBHE; i++)
    {
      double step=(eta_bins_HBHE[i+1]-eta_bins_HBHE[i])/totalMultiEta;
      for (int j=0; j<totalMultiEta; j++)
      {
        adjEtaBins[proj].push_back(eta_bins_HBHE[i]+step*j);
      }
    }
    adjEtaBins[proj].push_back(eta_bins_HBHE[eta_nbins_HBHE]);

    totalEtaBins[proj] = totalMultiEta*(eta_nbins_HBHE);
    totalPhiBins[proj] = granularityMultiPhi[proj] * granularityMultiECAL*HBHE_IPHI_NUM;

  }
	
  // Detector image switches
  doECALstitched = iConfig.getParameter<bool>("doECALstitched");
  doTracksAtECALstitchedPt = iConfig.getParameter<bool>("doTracksAtECALstitchedPt");
  doTracksAtECALadjPt = iConfig.getParameter<bool>("doTracksAtECALadjPt");
  doHBHEenergy = iConfig.getParameter<bool>("doHBHEenergy");
	
  //produces<float>("photonClasses").setBranchAlias("PhotonClass");
  produces<e2e::Frame1D>("EBenergy");
  produces<e2e::Frame1D>("HBHEenergyEB");
  produces<e2e::Frame1D>("TracksAtECALadj");
  produces<e2e::Frame1D>("TracksAtECALadjPtMax");
  produces<e2e::Frame3D>("DetFrames");
}

DetFrameProducer::~DetFrameProducer()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
DetFrameProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   //std::cout<<"New Event started"<<std::endl;
   using namespace edm;
   nTotal++;
   // ----- Apply event selection cuts ----- //
   
   //auto photon_classes = std::make_unique<float>(10.0);
   fillEB( iEvent, iSetup );
   std::unique_ptr<e2e::Frame1D> EBenergy_edm (new e2e::Frame1D(vEB_energy_));
	
   std::cout<<" >> Adding EB done "<<std::endl;
   std::cout<<" >> Size of EB Energy vector is: "<<std::move(EBenergy_edm).get()->size()<<std::endl;
   iEvent.put(std::move(EBenergy_edm),"EBenergy");
	
   e2e::Frame3D vDetFrames;
   e2e::Frame2D vHBHE_energy_reshaped (nDetFrameH, e2e::Frame1D (nDetFrameW,0.));
   e2e::Frame2D vECAL_energy_reshaped (nDetFrameH, e2e::Frame1D (nDetFrameW,0.));
   e2e::Frame2D vECAL_tracksPt_reshaped (nDetFrameH, e2e::Frame1D (nDetFrameW,0.));
   e2e::Frame2D vECALadj_tracksPt_reshaped (nDetFrameH, e2e::Frame1D (nDetFrameW,0.));
 
   if (doHBHEenergy){
   	fillHBHE (iEvent, iSetup );
	// reshape detector image arrays to 280x360
	for (unsigned int idx=0; idx<vHBHE_energy_.size(); idx++){
		vHBHE_energy_reshaped[int(idx/nDetFrameW)][idx%nDetFrameW]=vHBHE_energy_[idx];
	}
   	std::unique_ptr<e2e::Frame2D> HBHEenergy_edm (new e2e::Frame2D(vHBHE_energy_reshaped));
   	std::unique_ptr<e2e::Frame1D> HBHEenergyEB_edm (new e2e::Frame1D(vHBHE_energy_EB_));
   	std::cout<<" >> Size of HBHE Energy vector is: "<<std::move(HBHEenergy_edm).get()->size()<<std::endl;
   	std::cout<<" >> Size of EB HBHE Energy vector is: "<<std::move(HBHEenergyEB_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(HBHEenergy_edm),"HBHEenergy");
   	//iEvent.put(std::move(HBHEenergyEB_edm),"HBHEenergyEB");
	vDetFrames.push_back(vHBHE_energy_reshaped);
   }
   
   if (doECALstitched){
	fillECALstitched (iEvent, iSetup);
	// reshape detector image arrays to 280x360
	for (unsigned int idx=0; idx<vECAL_energy_.size(); idx++){
		vECAL_energy_reshaped[int(idx/nDetFrameW)][idx%nDetFrameW]=vECAL_energy_[idx];
	}
   	std::unique_ptr<e2e::Frame1D> ECALstitched_energy_edm (new e2e::Frame1D(vECAL_energy_));
   	std::cout<<" >> Size of Stitched ECAL Energy vector is: "<<std::move(ECALstitched_energy_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(ECALstitched_energy_edm), "ECALstitchedenergy");
   	vDetFrames.push_back(vECAL_energy_reshaped);
   }

   if (doTracksAtECALstitchedPt){
   	fillTracksAtECALstitched (iEvent, iSetup );
	// reshape detector image arrays to 280x360
	for (unsigned int idx=0; idx<vECAL_tracksPt_.size(); idx++){
		vECAL_tracksPt_reshaped[int(idx/nDetFrameW)][idx%nDetFrameW]=vECAL_tracksPt_[idx];
	}
   	std::unique_ptr<e2e::Frame1D> TracksECALstitchedPt_edm (new e2e::Frame1D(vECAL_tracksPt_));
   	std::cout<<" >> Size of Pt Tracks vector at Stitched ECAL is: "<<std::move(TracksECALstitchedPt_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(TracksECALstitchedPt_edm), "TracksAtECALstitchedPt");
   	vDetFrames.push_back(vECAL_tracksPt_reshaped);
   }

   if (doTracksAtECALadjPt){
   	for (unsigned int i=0;i<Nadjproj;i++)
   	{
     		fillTracksAtECALadjustable( iEvent, iSetup, i );
     		//fillTRKlayersAtECALadjustable( iEvent, iSetup, i );
   	}
	// reshape detector image arrays to 280x360
	for (unsigned int idx=0; idx<vECALadj_tracksPt_[0].size(); idx++){
		vECALadj_tracksPt_reshaped[int(idx/nDetFrameW)][idx%nDetFrameW]=vECALadj_tracksPt_[0][idx];
	}
   	std::cout<<" >> Number of TracksAtECALadjPt per event: "<<sizeof(vECALadj_tracksPt_)/sizeof(vECALadj_tracksPt_[0])<<std::endl;
   	std::cout<<" >> Number of TracksAtECALadj per event: "<<sizeof(vECALadj_tracksPt_)/sizeof(vECALadj_tracksPt_[0])<<std::endl;
   	std::cout<<" >> Number of TracksAtECALadjPtMax per event: "<<sizeof(vECALadj_tracksPt_max_)/sizeof(vECALadj_tracksPt_max_[0])<<std::endl;
   	std::cout<<" >> Sizes of TracksadjPt, Tracksadj and TracksadjPtMax are: "<<vECALadj_tracksPt_[0].size()<<", "<<vECALadj_tracks_[0].size()<<", "<<vECALadj_tracksPt_max_[0].size()<<std::endl;
   	std::unique_ptr<e2e::Frame1D> TracksECALadjPt_edm (new e2e::Frame1D(vECALadj_tracksPt_[0]));
   	std::cout<<" >> Size of Pt Tracks vector at ECAL adjustable is : "<<std::move(TracksECALadjPt_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(TracksECALadjPt_edm),"TracksAtECALadjPt");
   	std::unique_ptr<e2e::Frame1D> TracksECALadj_edm (new e2e::Frame1D(vECALadj_tracks_[0]));
   	std::cout<<" >> Size of Track vector at ECAL adjustable is : "<<std::move(TracksECALadj_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(TracksECALadj_edm),"TracksAtECALadj");
   	std::unique_ptr<e2e::Frame1D> TracksECALadjPt_max_edm (new e2e::Frame1D(vECALadj_tracksPt_max_[0]));
   	std::cout<<" >> Size of max Pt Track vector at ECAL adjustable is : "<<std::move(TracksECALadjPt_max_edm).get()->size()<<std::endl;
   	//iEvent.put(std::move(TracksECALadjPt_max_edm),"TracksAtECALadjPtMax");
   	vDetFrames.push_back(vECALadj_tracksPt_reshaped);
   }	
   std::unique_ptr<e2e::Frame3D> vDetFrames_edm (new e2e::Frame3D(vDetFrames));
   iEvent.put(std::move(vDetFrames_edm), "DetFrames");
   
   std::cout<<" >> Added EB, HBHE, HBHE_EB, ECALstitched, TracksAtECALstitchedPt and TracksAtECALadjPt to edm root file"<<std::endl;
   
   nPassed++;
   std::cout<<std::endl;
   return;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
DetFrameProducer::beginStream(edm::StreamID)
{
 nTotal = 0;
 nPassed = 0;
 std::cout<<"'DetFrameProducer' Stream began"<<std::endl;
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
DetFrameProducer::endStream() {
 std::cout << "'DetFrameProducer' selected: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void
ProducerTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
ProducerTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
ProducerTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
ProducerTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DetFrameProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(DetFrameProducer);
