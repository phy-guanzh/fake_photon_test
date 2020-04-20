// -*- C++ -*-
//
// Package:    photonAna/photonAnalyzer
// Class:      photonAnalyzer
//
/**\class photonAnalyzer photonAnalyzer.cc photonAna/photonAnalyzer/plugins/photonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Thu, 28 Nov 2019 03:34:32 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class photonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit photonAnalyzer(const edm::ParameterSet&);
      ~photonAnalyzer();

//      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      void setDummyValues();


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
	float EAch(float x); 
	float EAnh(float x);
	float EApho(float x);

      // ----------member data ---------------------------
//      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<pat::Photon>> photonToken_;  
	edm::EDGetTokenT<double> rhoToken_;
	edm::Handle< double >  rho_;
	Double_t photon_pt;
	Double_t photon_eta;
	Double_t photon_phi;
	Double_t photon_e;
	Double_t HoverE;
	Double_t sieie;
	Double_t chIso;
	Double_t neuIso;
	Double_t phoIso;
	Bool_t passEleVeto;
	TTree* outTree_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
photonAnalyzer::photonAnalyzer(const edm::ParameterSet& iConfig)
 :
//  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
  photonToken_(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("Photons")))
{
   //now do what ever initialization is needed
  rhoToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
   edm::Service<TFileService> fs;
   outTree_ = fs->make<TTree>("PhotonCandidates","Photon Candidates");
   outTree_->Branch("photon_pt", &photon_pt, "photon_pt/D");
   outTree_->Branch("photon_eta", &photon_eta, "photon_eta/D");
   outTree_->Branch("photon_phi", &photon_phi, "photon_phi/D");
   outTree_->Branch("photon_e", &photon_e, "photon_e/D");
   outTree_->Branch("HoverE", &HoverE, "HoverE/D");
   outTree_->Branch("sieie", &sieie, "sieie/D");
   outTree_->Branch("chIso", &chIso, "chIso/D");
   outTree_->Branch("neuIso", &neuIso, "neuIso/D");
   outTree_->Branch("phoIso", &phoIso, "phoIso/D");
   outTree_->Branch("passEleVeto", &passEleVeto, "passEleVeto/O");

}


photonAnalyzer::~photonAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
photonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 //   Handle<TrackCollection> tracks;
 //   iEvent.getByToken(tracksToken_, tracks);
 //   for(TrackCollection::const_iterator itTrack = tracks->begin();
 //       itTrack != tracks->end();
 //       ++itTrack) {
 //     // do something with track parameters, e.g, plot the charge.
 //     // int charge = itTrack->charge();
 //   }
	setDummyValues();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
   Handle<edm::View<pat::Photon>> photons;
   iEvent.getByToken(photonToken_, photons);

   iEvent.getByToken(rhoToken_ ,rho_ );
   double rhoVal_;
   rhoVal_=-99.;
   rhoVal_ = *rho_;
   // The reco-level photon candidates have been sorted by pT, we just select the leading-pT one in a single event for simplicity 
   if(photons->size()<1) {outTree_->Fill();return;}
   const auto pho = photons->ptrAt(0);
   photon_pt = pho->pt();
   photon_eta = pho->superCluster()->eta();
   photon_phi = pho->superCluster()->phi();
   photon_e = pho->energy();
   HoverE = pho->hadTowOverEm();
   sieie = pho->full5x5_sigmaIetaIeta();
   // PFIso_corrected = max(PFIso - rho*EA, 0.)
   chIso = pho->chargedHadronIso();
   chIso = std::max(0.0, chIso - rhoVal_*EAch(fabs(photon_eta)));
   neuIso = pho->neutralHadronIso();
   neuIso = std::max(0.0, neuIso - rhoVal_*EAnh(fabs(photon_eta)));
   phoIso = pho->photonIso();
   phoIso = std::max(0.0, phoIso - rhoVal_*EApho(fabs(photon_eta)));
   passEleVeto = pho->passElectronVeto();
   outTree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void
photonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
photonAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void
//photonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//  //The following says we do not know what parameters are allowed so do no validation
//  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
//
//  //Specify that only 'tracks' is allowed
//  //To use, remove the default given above and uncomment below
//  //ParameterSetDescription desc;
//  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
//  //descriptions.addDefault(desc);
//}

void
photonAnalyzer::setDummyValues(){
	photon_pt=-99.;
        photon_eta=-99.;
        photon_phi=-99.;
        photon_e=-99.;
        HoverE=-99.;
        sieie=-99.;
        chIso=-99.;
        neuIso=-99.;
        phoIso=-99.;
        passEleVeto=-99;
}

float photonAnalyzer::EAch( float x){
	float EA = 0.0112;
	if(x>1.0)   EA = 0.0108;
	if(x>1.479) EA = 0.0106;
	if(x>2.0)   EA = 0.01002;
	if(x>2.2)   EA = 0.0098;
	if(x>2.3)   EA = 0.0089;
	if(x>2.4)   EA = 0.0087;
	return EA;
}
float photonAnalyzer::EAnh( float x){
	float EA = 0.0668;
	if(x>1.0)   EA = 0.1054;
	if(x>1.479) EA = 0.0786;
	if(x>2.0)   EA = 0.0233;
	if(x>2.2)   EA = 0.0078;
	if(x>2.3)   EA = 0.0028;
	if(x>2.4)   EA = 0.0137;
	return EA;
}
float photonAnalyzer::EApho( float x){
	float EA = 0.1113;
	if(x>1.0)   EA = 0.0953;
	if(x>1.479) EA = 0.0619;
	if(x>2.0)   EA = 0.0837;
	if(x>2.2)   EA = 0.1070;
	if(x>2.3)   EA = 0.1212;
	if(x>2.4)   EA = 0.1466;
	return EA;
}

//define this as a plug-in
DEFINE_FWK_MODULE(photonAnalyzer);
