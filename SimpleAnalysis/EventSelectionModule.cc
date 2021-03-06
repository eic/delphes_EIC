#include "EventSelectionModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"
#include "TreeHandler.h"

#include <iostream>
#include <iomanip>  
#include <fstream>

EventSelectionModule::EventSelectionModule(ExRootTreeReader* data)
  : Module(data)
{

}

EventSelectionModule::~EventSelectionModule()
{

}

void EventSelectionModule::initialize()
{
  TreeHandler *tree_handler = tree_handler->getInstance();

  if (tree_handler->getTree() != nullptr) {

    _jet_n = 0;
    _jet_pt = std::vector<float>();
    _jet_eta = std::vector<float>();
    _jet_flavor = std::vector<int>();
    _jet_sip3dtag = std::vector<int>();
    _jet_ktag = std::vector<int>();
    _jet_etag = std::vector<int>();
    _jet_mutag = std::vector<int>();

    tree_handler->getTree()->Branch("jet_n",        &_jet_n, "jet_n/I");
    tree_handler->getTree()->Branch("jet_pt",       "std::vector<float>", &_jet_pt);
    tree_handler->getTree()->Branch("jet_eta",      "std::vector<float>", &_jet_eta);
    tree_handler->getTree()->Branch("jet_flavor",   "std::vector<int>", &_jet_flavor);
    tree_handler->getTree()->Branch("jet_sip3dtag", "std::vector<int>", &_jet_sip3dtag);
    tree_handler->getTree()->Branch("jet_ktag",     "std::vector<int>", &_jet_ktag);
    tree_handler->getTree()->Branch("jet_etag",     "std::vector<int>", &_jet_etag);
    tree_handler->getTree()->Branch("jet_mutag",    "std::vector<int>", &_jet_mutag);

    _jet_Ks_mass = std::vector<float>();
    _jet_Ks_p = std::vector<float>();
    _jet_Ks_flightlength = std::vector<float>();
    _jet_Ks_sumpt = std::vector<float>();
    _jet_K_sumpt = std::vector<float>();

    tree_handler->getTree()->Branch("jet_Ks_mass",          "std::vector<float>", &_jet_Ks_mass);
    tree_handler->getTree()->Branch("jet_Ks_p",             "std::vector<float>", &_jet_Ks_p);
    tree_handler->getTree()->Branch("jet_Ks_flightlength",  "std::vector<float>", &_jet_Ks_flightlength);
    tree_handler->getTree()->Branch("jet_Ks_sumpt",         "std::vector<float>", &_jet_Ks_sumpt);
    tree_handler->getTree()->Branch("jet_K_sumpt",          "std::vector<float>", &_jet_K_sumpt);


    _charmjet_n = 0;
    _charmjet_pt = std::vector<float>();
    _charmjet_eta = std::vector<float>();
    tree_handler->getTree()->Branch("charmjet_n", &_charmjet_n, "charmjet_n/I");
    tree_handler->getTree()->Branch("charmjet_pt", "std::vector<float>", &_charmjet_pt);
    tree_handler->getTree()->Branch("charmjet_eta", "std::vector<float>", &_charmjet_eta);

    _met_et = 0.0;
    tree_handler->getTree()->Branch("met_et", &_met_et, "met_et/F");

    _bjorken_x = 0.0;
    _bjorken_Q2 = 0.0;
    _bjorken_y = 0.0;
    _jb_x = 0.0;
    _jb_Q2 = 0.0;
    tree_handler->getTree()->Branch("bjorken_x", &_bjorken_x, "bjorken_x/F");
    tree_handler->getTree()->Branch("bjorken_Q2", &_bjorken_Q2, "bjorken_Q2/F");
    tree_handler->getTree()->Branch("bjorken_y", &_bjorken_y, "bjorken_y/F");
    tree_handler->getTree()->Branch("jb_x", &_jb_x, "jb_x/F");
    tree_handler->getTree()->Branch("jb_Q2", &_jb_Q2, "jb_Q2/F");

  }


  // Initialize the cut flow

  _cut_flow["1: All events"] = 0;
  _cut_flow["2: MET > 10 GeV"] = 0;
  _cut_flow["3: Fiducial Jets >= 1"] = 0;
  _cut_flow["4: Charm Jet == 1"] = 0;

  // Global variables
  _mpi = TDatabasePDG().GetParticle(211)->Mass();



}

void EventSelectionModule::finalize()
{
  ofstream csvfile;
  csvfile.open("cut_flow.csv");
  
  csvfile << "Cut,Yield" << std::endl;

  for (auto& [cut,yield] : _cut_flow) {
    csvfile << "\"" << cut << "\"," << int(yield) << std::endl;
  }

  csvfile.close();

}


bool EventSelectionModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();

  // Compute global DIS variables
  auto dis_variables = DISVariables(getGenParticles());
  _bjorken_x = dis_variables["x"];
  _bjorken_Q2 = dis_variables["Q2"];
  _bjorken_y = dis_variables["y"];

  auto jb_variables = DISJacquetBlondel(getEFlowTracks(), getElectrons(), getPhotons(), getNeutralHadrons());
  _jb_x = jb_variables["x_JB"];
  _jb_Q2 = jb_variables["Q2_JB"];


  // Initialize output variables
  // _charmjet_pt = std::vector<float>();
  // _charmjet_eta = std::vector<float>();
  _charmjet_pt.clear();
  _charmjet_eta.clear();
  _charmjet_n = _charmjet_pt.size();

  // Cut flow
  _cut_flow["1: All events"] += 1;
  bool passed = true;

  // Get the MET object and use it
  MissingET* MET = nullptr;
  for (int imet = 0; imet < getMET()->GetEntries(); imet++) {
    MET = static_cast<MissingET*>(getMET()->At(imet));
  }
  
  if (MET == nullptr) {
    passed = false;
  }

  _met_et = MET->MET;
  
  if (passed == true && MET->MET > 10.0) {
    _cut_flow["2: MET > 10 GeV"] += 1;
  } else {
    passed = false;
  }


  // If event contains at least 1 jet


  _jet_n = 0;
  // _jet_pt = std::vector<float>();
  // _jet_eta = std::vector<float>();
  // _jet_flavor = std::vector<int>();
  // _jet_sip3dtag = std::vector<int>();
  // _jet_ktag = std::vector<int>();
  // _jet_etag = std::vector<int>();
  // _jet_mutag = std::vector<int>();
  _jet_pt.clear();
  _jet_eta.clear();
  _jet_flavor.clear();
  _jet_sip3dtag.clear();
  _jet_ktag.clear();
  _jet_etag.clear();
  _jet_mutag.clear();

  _jet_Ks_mass.clear();
  _jet_Ks_p.clear();
  _jet_Ks_flightlength.clear();
  _jet_Ks_sumpt.clear();
  _jet_K_sumpt.clear();
  
  auto tracks = getTracks();


  bool use_kaons = false;
  if (DataStore->find("Kaons") != DataStore->end()) {
    // store the number of kaons in the jets
    use_kaons = true;
  }

  bool use_electrons = false;
  if (DataStore->find("Electrons") != DataStore->end()) {
    // store the number of electrons in the jets
    use_electrons = true;
  }

  bool use_muons = false;
  if (DataStore->find("Muons") != DataStore->end()) {
    // store the number of muons in the jets
    use_muons = true;
  }


  // Loop over tracks in the event; make opposite-sign pairs; compute mass and call them Ks if within some window
  // of the Ks0 mass.

  // auto particles = getGenParticles();

  std::vector<Track*> all_tracks;
  //std::vector<GenParticle*> all_tracks;
  //std::vector<TLorentzVector> all_tracks;
  for (auto obj_track : *tracks) 
    {
      auto track = static_cast<Track*>(obj_track);
      if (track->PT < 0.1)
	continue;
      all_tracks.push_back( track );
    }


  // Build all Ks candidates
  std::vector<Candidate> Ks_candidates;

  for (int i = 0; i < all_tracks.size(); i++)
    {
      for (int j = i + 1; j < all_tracks.size(); j++) 
	{
	  auto track1 = all_tracks[i];
	  auto track2 = all_tracks[j];

	  if (track1->Charge * track2->Charge != -1)
	     continue;

	  // Treat the tracks under the charged pion hypothesis
	  auto track1P4 = TLorentzVector();
	  track1P4.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, _mpi);
	  auto track2P4 = TLorentzVector();
	  track2P4.SetPtEtaPhiM(track2->PT, track2->Eta, track2->Phi, _mpi);
	      
	  // std::cout << "===========================================================" << std::endl;
	  // track1->P4().Print();
	  // track2->P4().Print();
	  
	  //TLorentzVector Ks_candidate = track1->P4() + track2->P4();
	  TLorentzVector Ks_candidate = track1P4 + track2P4;
	      
	  if ( TMath::Abs(Ks_candidate.M() - 0.497) > 0.250 )
	    continue;

	  // Spatial coincidence requirement
	  TVector3 track1_POCA(track1->X, track1->Y, track1->Z);
	  TVector3 track2_POCA(track2->X, track2->Y, track2->Z);
	  
	  auto intertrack_displacement = track1_POCA - track2_POCA;
	  float intertrack_distance = intertrack_displacement.Mag();

	  float track1_d0err = track1->ErrorD0;
	  float track1_z0err = track1->ErrorDZ;
	  float track2_d0err = track2->ErrorD0;
	  float track2_z0err = track2->ErrorDZ;

	  float err_3D = TMath::Sqrt( TMath::Power(track1_d0err,2.0) + TMath::Power(track1_z0err,2.0) +
				      TMath::Power(track2_d0err,2.0) + TMath::Power(track2_z0err,2.0) );
	  
	  float intertrack_distance_signif = intertrack_distance/err_3D;

	  // std::cout << " Intertrack distance, error, significance: " << intertrack_distance << " mm" 
	  // 	    << ", " << err_3D 
	  // 	    << ", " << intertrack_distance_signif << std::endl;

	  // To be from a common decay, their displacement significance should be small

	  if (intertrack_distance_signif > 1.5)
	    continue;

	  // build a new Candidate
	  Candidate Ks;
	  Ks.PID = 310;
	  Ks.Mass = Ks_candidate.M();
	  Ks.Momentum = Ks_candidate;
	  Ks.Position = TLorentzVector(track2_POCA + 0.5*intertrack_displacement, 0.0); // midpoint between POCA of two tracks
	  

	      
	  Ks_candidates.push_back(Ks);
	}
    }

  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      // Take first jet
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);

      _jet_pt.push_back( jet->PT );
      _jet_eta.push_back( jet->Eta );
      _jet_flavor.push_back( jet->Flavor );
      _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *tracks, 3.75, 1.00, 2.0) );
      if (use_electrons) {
	_jet_etag.push_back( Tagged_Electron(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Electrons"]), 3.0, 1.0, 1) );
      } else { 
	_jet_etag.push_back(0.0);
      }


      if (use_muons) {
	_jet_mutag.push_back( Tagged_Muon(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Muons"]), 3.0, 1.0, 1) );
      } else {
	_jet_mutag.push_back( 0.0 );
      }

      if (use_kaons) {
	_jet_ktag.push_back( Tagged_Kaon(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Kaons"]), 3.0, 1.0, 1) );
      } else {
	_jet_ktag.push_back( 0.0 );
      }

      Candidate best_Ks;
      TVector3 Ks_sumpt;
      for (auto Ks_candidate : Ks_candidates) {
	if (Ks_candidate.Position.Rho() < 5) // 5mm minimum displacement from IP 
	  continue;
	if (Ks_candidate.Momentum.DeltaR( jet->P4() ) < 0.5) {
	  
	  Ks_sumpt += Ks_candidate.Momentum.Vect();

	  if (_jet_Ks_mass.size() < ijet+1) {
	    _jet_Ks_mass.push_back( Ks_candidate.Mass );
	    _jet_Ks_p.push_back( Ks_candidate.Momentum.Rho() );
	    _jet_Ks_flightlength.push_back( Ks_candidate.Position.Rho() );
	    best_Ks = Ks_candidate;
	  } else {
	    if (Ks_candidate.Position.Rho() > best_Ks.Position.Rho()) {
	      _jet_Ks_mass[ijet] = Ks_candidate.Mass;
	      _jet_Ks_p[ijet] = Ks_candidate.Momentum.Rho();
	      _jet_Ks_flightlength[ijet] = Ks_candidate.Position.Rho();
	      best_Ks = Ks_candidate;
	    }
	  }
	}
      }
      _jet_Ks_sumpt.push_back( Ks_sumpt.Perp() );

      // handle charged kaons
      TVector3 K_sumpt;
      if (use_kaons) {
	auto kaon_candidates = std::any_cast<std::vector<Track*>>((*DataStore)["Kaons"]);
	for (auto kaon : kaon_candidates) {
	  if (kaon->P4().DeltaR(jet->P4()) < 0.5) {
	    K_sumpt += kaon->P4().Vect();
	  }
	}
      }
      _jet_K_sumpt.push_back( K_sumpt.Perp() );


    }

  _jet_n = _jet_pt.size();

  std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });

  if (passed == true && fiducial_jets.size() > 0) {
    _cut_flow["3: Fiducial Jets >= 1"] += 1;
  } else {
    passed = false;
  }

  //std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor == 4); });
  
  std::vector<Jet*> charmJets;
  if (DataStore->find("CharmJets") != DataStore->end()) {
    charmJets = std::any_cast<std::vector<Jet*>>((*DataStore)["CharmJets"]);
  }

  if (passed == true && charmJets.size() > 0) {
    _cut_flow["4: Charm Jet == 1"] += 1;
  } else {
    passed = false;
  }


  // Store charm jet information
  if (passed) {
    for (auto jet : charmJets) {
      _charmjet_pt.push_back( jet->PT );
      _charmjet_eta.push_back( jet->Eta );
    }
    _charmjet_n = _charmjet_pt.size();
  }


  return true;
}




