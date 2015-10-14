// -*- C++ -*-
//
// input: L1 TkTracks and  L1MuonParticleExtended (standalone with component details)
// match the two and produce a collection of L1TkMuonParticle
// eventually, this should be made modular and allow to swap out different algorithms


#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkMuonParticleFwd.h"


// for L1Tracks:
#include "SimDataFormats/SLHC/interface/StackedTrackerTypes.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
 #include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMath.h"

// system include files
#include <memory>
#include <string>
#include <TTree.h>

using namespace l1extra ;

//
// class declaration
//

class L1TkMuonFromExtendedProducer : public edm::EDProducer {
public:
  
  typedef TTTrack< Ref_PixelDigi_ >                     L1TkTrackType;
  typedef std::vector< L1TkTrackType >        L1TkTrackCollectionType;
  
  struct PropState { //something simple, imagine it's hardware emulation
    PropState() : 
      pt(-99),  eta(-99), phi(-99),
      sigmaPt(-99),  sigmaEta(-99), sigmaPhi(-99),
      valid(false) {}
    float pt;
    float eta;
    float phi;
    float sigmaPt;
    float sigmaEta;
    float sigmaPhi;
    bool valid;

  };

  explicit L1TkMuonFromExtendedProducer(const edm::ParameterSet&);
  ~L1TkMuonFromExtendedProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  PropState propagateToGMT(const L1TkTrackType& l1tk) const;

  //configuration (preserving structure of L1TkMuonFromExtendedProducer
  edm::InputTag L1MuonsInputTag_;
  edm::InputTag L1TrackInputTag_;	 
  
  float ETAMIN_;
  float ETAMAX_;
  float ZMAX_;             // |z_track| < ZMAX in cm
  float CHI2MAX_;
  float PTMINTRA_;
  //  float DRmax_;
  int nStubsmin_ ;         // minimum number of stubs   
  //  bool closest_ ;
  bool correctGMTPropForTkZ_;

  bool use5ParameterFit_;

  float tp_pt, tp_eta, tp_phi, tp_z;
  float tk_pt, tk_eta, tk_phi, tk_z;
  float mu_pt, mu_eta, mu_phi, mu_z;
  int Genuine, Combinatoric;
  int nTracks;
  TTree*Amtracks;


} ;


//
// constructors and destructor
//
L1TkMuonFromExtendedProducer::L1TkMuonFromExtendedProducer(const edm::ParameterSet& iConfig)
{


   L1MuonsInputTag_ = iConfig.getParameter<edm::InputTag>("L1MuonsInputTag");
   L1TrackInputTag_ = iConfig.getParameter<edm::InputTag>("L1TrackInputTag");

   ETAMIN_ = (float)iConfig.getParameter<double>("ETAMIN");
   ETAMAX_ = (float)iConfig.getParameter<double>("ETAMAX");
   ZMAX_ = (float)iConfig.getParameter<double>("ZMAX");
   CHI2MAX_ = (float)iConfig.getParameter<double>("CHI2MAX");
   PTMINTRA_ = (float)iConfig.getParameter<double>("PTMINTRA");
   //   DRmax_ = (float)iConfig.getParameter<double>("DRmax");
   nStubsmin_ = iConfig.getParameter<int>("nStubsmin");
   //   closest_ = iConfig.getParameter<bool>("closest");

   correctGMTPropForTkZ_ = iConfig.getParameter<bool>("correctGMTPropForTkZ");

   use5ParameterFit_     = iConfig.getParameter<bool>("use5ParameterFit");
   produces<L1TkMuonParticleCollection>();
}

L1TkMuonFromExtendedProducer::~L1TkMuonFromExtendedProducer() {
}
void L1TkMuonFromExtendedProducer::beginJob(){
edm::Service<TFileService> fs;
    Amtracks = fs->make<TTree>("Amtracks", "");
    Amtracks->Branch("tp_eta", &tp_eta, "tp_eta/F");
    Amtracks->Branch("tp_phi", &tp_phi, "tp_phi/F");
    Amtracks->Branch("tp_pt", &tp_pt, "tp_pt/F");
    Amtracks->Branch("tp_z", &tp_z, "tp_z/F");
    Amtracks->Branch("tk_eta", &tk_eta, "tk_eta/F");
    Amtracks->Branch("tk_phi", &tk_phi, "tk_phi/F");
    Amtracks->Branch("tk_pt", &tk_pt, "tk_pt/F");
    Amtracks->Branch("tk_z", &tk_z, "tk_z/F");
    Amtracks->Branch("Genuine", &Genuine, "Genuine/I");
 Amtracks->Branch("Combinatoric", &Combinatoric, "Combinatoric/I");
    Amtracks->Branch("mu_eta", &mu_eta, "mu_eta/F");
    Amtracks->Branch("mu_phi", &mu_phi, "mu_phi/F");
    Amtracks->Branch("mu_pt", &mu_pt, "mu_pt/F");
    Amtracks->Branch("mu_z", &mu_z, "mu_z/F");
}
// ------------ method called to produce the data  ------------
void
L1TkMuonFromExtendedProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
 edm::Handle< TTStubAssociationMap< Ref_PixelDigi_ > > MCTruthTTStubHandle;
 iEvent.getByLabel("TTStubAssociatorFromPixelDigis", "StubAccepted", MCTruthTTStubHandle);

 edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
iEvent.getByLabel("mix", "MergedTrackTruth", TrackingParticleHandle);

std::vector< TrackingParticle >::const_iterator iterTP;
tp_eta=999;
tp_phi=999;
tp_pt=-1;
tp_z=999;
 for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    if(iterTP->pdgId()!=13)continue;
    tp_eta=iterTP->eta();
    tp_phi=iterTP->phi();
    tp_pt=iterTP->pt();
    tp_z=iterTP->vz();
}
float dRMin=99.;
tk_eta=999;
tk_phi=999;
tk_pt=-1;
tk_z=999;
mu_eta=999;
mu_phi=999;
mu_pt=-1;
mu_z=999; 
  
  edm::Handle<L1TkTrackCollectionType> l1tksH;
  iEvent.getByLabel(L1TrackInputTag_, l1tksH);
  const L1TkTrackCollectionType& l1tks = *l1tksH.product();
        Genuine=0;
        Combinatoric=0;

for (auto l1tk : l1tks ){

      float deta= fabs(l1tk.getMomentum(4).eta()-tp_eta);
      float dphi= acos(cos(l1tk.getMomentum(4).phi()-tp_phi));
      float dR=sqrt((deta*deta)+(dphi*dphi));
      if(dR<dRMin){
	dRMin=dR;
	tk_eta=l1tk.getMomentum(4).eta();
	tk_phi=l1tk.getMomentum(4).phi();
	tk_pt=l1tk.getMomentum(4).perp();
        tk_z=l1tk.getPOCA(4).z();
      int l1tk_nstubs = l1tk.getStubRefs().size();
	Genuine=0;
        Combinatoric=0;	
        for(int s=0; s<l1tk_nstubs; ++s){
                if(MCTruthTTStubHandle.isValid()){
                          bool isGenuine = MCTruthTTStubHandle->isGenuine(l1tk.getStubRefs()[s]);
			  bool isComb = MCTruthTTStubHandle->isCombinatoric(l1tk.getStubRefs()[s]);
                                      if(isGenuine)++Genuine;
				      if(isComb)++Combinatoric;
                                              }
                                                  }
      }
}
 
//if(MCTruthTTStubHandle.isValid())


  std::auto_ptr<L1TkMuonParticleCollection> tkMuons(new L1TkMuonParticleCollection);
  
  edm::Handle<L1MuonParticleExtendedCollection> l1musH;
  iEvent.getByLabel(L1MuonsInputTag_, l1musH);
  const L1MuonParticleExtendedCollection& l1mus = *l1musH.product(); 
  
  L1TkMuonParticleCollection l1tkmuCands;
  l1tkmuCands.reserve(l1mus.size()*4); //can do more if really needed
  nTracks=l1tks.size();
  std::cout<<"Track Collection Size "<<l1tks.size()<<std::endl;
  std::cout<<"Extended Collection Size "<<l1mus.size()<<std::endl;
  int imu = -1;
  for (auto l1mu : l1mus){
    imu++;
    L1MuonParticleExtendedRef l1muRef(l1musH, imu);

    float l1mu_eta = l1mu.eta();
    float l1mu_phi = l1mu.phi();
    
    float l1mu_feta = fabs( l1mu_eta );
    if (l1mu_feta < ETAMIN_) continue;
    if (l1mu_feta > ETAMAX_) continue;

    // can skip quality cuts at the moment, keep bx=0 req
    if (l1mu.bx() != 0) continue;

    unsigned int l1mu_quality = l1mu.quality();

    const auto& gmtCand = l1mu.gmtMuonCand();
    if (!gmtCand.empty()){
      //some selections here
      //currently the input can be a merge from different track finders
      //so, the GMT cand may be missing
    }

    const auto& dtCand  = l1mu.dtCand();
    if (!dtCand.empty()){
      // something can be called from here
    }

    const auto& cscCand = l1mu.cscCand();
    if (!cscCand.empty()){
      //apply something specific here
    }

    const auto& rpcCand = l1mu.rpcCand();
    if (!rpcCand.empty()){
      //apply something specific here
    }

    float drmin = 999;
    float ptmax = -1;
    if (ptmax < 0) ptmax = -1;	// dummy
    
    PropState matchProp;
    int match_idx = -1;
    int il1tk = -1;

    std::cout<<"have a gmt, look for a match "<<std::endl;
    for (auto l1tk : l1tks ){
      il1tk++;

      unsigned int nPars = 4;
      if (use5ParameterFit_) nPars = 5;
      float l1tk_pt = l1tk.getMomentum(nPars).perp();
      // std::cout<<"Pt "<<l1tk_pt<<std::endl;
      if (l1tk_pt < PTMINTRA_) continue;
      //std::cout<<"Pt Filter "<<std::endl;
      float l1tk_z  = l1tk.getPOCA(nPars).z();
      if (fabs(l1tk_z) > ZMAX_) continue;
      //std::cout<<"Z Filter "<<std::endl;
      //float l1tk_chi2 = l1tk.getChi2(nPars);
     // std::cout<<"Chi2 "<<l1tk.getChi2(nPars)<<std::endl;
      //if (l1tk_chi2 > CHI2MAX_) continue;
      //std::cout<<"chi2 Filter "<<std::endl;
      int l1tk_nstubs = l1tk.getStubRefs().size();
//std::vector< edm::Ref< TTStub > >SRefs=l1tk.getStubRefs();
/*
	for(int s=0; s<l1tk_nstubs; ++s){
	if(MCTruthTTStubHandle.isValid()){
	   //bool isGenuine = MCTruthTTStubHandle->isGenuine(l1tk.getStubRefs()[s]);
	 //  if(isGenuine)std::cout<<"Genuine Stub in Track "<<std::endl;
	}
    }
*/
    if ( l1tk_nstubs < nStubsmin_) continue;

      float l1tk_eta = l1tk.getMomentum(nPars).eta();
      float l1tk_phi = l1tk.getMomentum(nPars).phi();


      float dr2 = deltaR2(l1mu_eta, l1mu_phi, l1tk_eta, l1tk_phi);
     // std::cout<<"Loop over Tracks dRValue "<<dr2<<std::endl;
      if (dr2 > 0.3) continue;

      PropState pstate = propagateToGMT(l1tk);
      
      if (!pstate.valid) continue;
    //  else       std::cout<<"Propagator is Valid "<<std::endl;
      float dr2prop = deltaR2(l1mu_eta, l1mu_phi, pstate.eta, pstate.phi);
      std::cout<<"Match Radius "<<dr2prop<<std::endl;
      if (dr2prop < drmin){
	drmin = dr2prop;
	match_idx = il1tk;
	matchProp = pstate;
      }
    }// over l1tks
    
    LogDebug("MYDEBUG")<<"matching index is "<<match_idx;
    std::cout<<"matching index is "<<match_idx<<std::endl;
    if (match_idx >= 0){
      const L1TkTrackType& matchTk = l1tks[match_idx];
        
      float etaCut = 3.*sqrt(l1mu.sigmaEta()*l1mu.sigmaEta() + matchProp.sigmaEta*matchProp.sigmaEta);
      float phiCut = 4.*sqrt(l1mu.sigmaPhi()*l1mu.sigmaPhi() + matchProp.sigmaPhi*matchProp.sigmaPhi);

      float dEta = std::abs(matchProp.eta - l1mu.eta());
      float dPhi = std::abs(deltaPhi(matchProp.phi, l1mu.phi()));

      //LogDebug("MYDEBUG")
	//std::cout<<"match details: prop "<<matchProp.pt<<" "<<matchProp.eta<<" "<<matchProp.phi
	//		 <<" mutk "<<l1mu.pt()<<" "<<l1mu.eta()<<" "<<l1mu.phi()<<" delta "<<dEta<<" "<<dPhi<<" cut "<<etaCut<<" "<<phiCut<<std::endl;
      if (dEta < etaCut && dPhi < phiCut){
	Ptr< L1TkTrackType > l1tkPtr(l1tksH, match_idx);

	unsigned int nPars = 4;
	if (use5ParameterFit_) nPars = 5;
	auto p3 = matchTk.getMomentum(nPars);
	float p4e = sqrt(0.105658369*0.105658369 + p3.mag2() );

	math::XYZTLorentzVector l1tkp4(p3.x(), p3.y(), p3.z(), p4e);

	auto tkv3=matchTk.getPOCA(nPars);
	math::XYZPoint v3(tkv3.x(), tkv3.y(), tkv3.z());
	float trkisol = -999;
	int l1tk_q = matchTk.getRInv(nPars)>0? 1: -1;
        mu_eta=l1tkp4.eta();
	mu_phi=l1tkp4.phi();
	mu_pt=l1tkp4.pt();
	mu_z=tkv3.z();
	L1TkMuonParticle l1tkmu(reco::LeafCandidate(l1tk_q, l1tkp4, v3, -13*l1tk_q ));
	l1tkmu.setTrkPtr(l1tkPtr);
	l1tkmu.setMuExtendedRef(l1muRef);
	l1tkmu.setQuality(l1mu_quality);
	l1tkmu.setTrkIsol(trkisol);

	// EP: add the zvtx information
	l1tkmu.setTrkzVtx( (float)tkv3.z() );

	tkMuons->push_back(l1tkmu);
      }
    }

  }//over l1mus
  
  Amtracks->Fill(); 
  iEvent.put( tkMuons );
  
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TkMuonFromExtendedProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


L1TkMuonFromExtendedProducer::PropState L1TkMuonFromExtendedProducer::propagateToGMT(const L1TkMuonFromExtendedProducer::L1TkTrackType& tk) const {
  auto p3 = tk.getMomentum();
  float tk_pt = p3.perp();
  float tk_p = p3.mag();
  float tk_eta = p3.eta();
  float tk_aeta = std::abs(tk_eta);
  float tk_phi = p3.phi();
  float tk_q = tk.getRInv()>0? 1.: -1.;
  float tk_z  = tk.getPOCA().z();
  if (!correctGMTPropForTkZ_) tk_z = 0;
  else std::cout<<"Z Correction "<<std::endl;
    std::cout<<"P "<<tk_p<<" PT "<<tk_pt<<" eta "<<tk_eta<<std::endl;
    std::cout<<"rinv "<<tk.getRInv()<<" Phi "<<tk_phi<<std::endl;
  L1TkMuonFromExtendedProducer::PropState dest;
  if (tk_p<3.5 ) return dest;
  if (tk_aeta <1.1 && tk_pt < 3.5) return dest;
  if (tk_aeta > 2.5) return dest;

  //0th order:
  dest.valid = true;

  float dzCorrPhi = 1.;
  float deta = 0;
  float etaProp = tk_aeta;

  if (tk_aeta < 1.1){
    etaProp = 1.1;
    deta = tk_z/550./cosh(tk_aeta);
  } else {
    float delta = tk_z/850.; //roughly scales as distance to 2nd station
    if (tk_eta > 0) delta *=-1;
    dzCorrPhi = 1. + delta;

    float zOzs = tk_z/850.;
    if (tk_eta > 0) deta = zOzs/(1. - zOzs);
    else deta = zOzs/(1.+zOzs);
    deta = deta*tanh(tk_eta);
  }
  float resPhi = tk_phi - 1.464*tk_q*cosh(1.7)/cosh(etaProp)/tk_pt*dzCorrPhi - M_PI/144.;
  if (resPhi > M_PI) resPhi -= 2.*M_PI;
  if (resPhi < -M_PI) resPhi += 2.*M_PI;

  dest.eta = tk_eta + deta;
  dest.phi = resPhi;
  dest.pt = tk_pt; //not corrected for eloss

  dest.sigmaEta = 0.100/tk_pt; //multiple scattering term
  dest.sigmaPhi = 0.106/tk_pt; //need a better estimate for these
  return dest;
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1TkMuonFromExtendedProducer);



