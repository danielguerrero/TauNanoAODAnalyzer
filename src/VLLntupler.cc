#include "VLLntupler.h"
#include "EventAnalyzer.h"
#include "JetTree.h"
#include <stdlib.h> 
#include "JetCorrectionUncertainty.h"
#include "JetCorrectorParameters.h"
//C++ includes
//ROOT includes
#include "TH1F.h"
#include "TFile.h"
#include <chrono>

//using namespace correction;
using namespace std;

const float ELE_MASS = 0.000511;
const float MU_MASS  = 0.105658;

//struct greater_than_pt
//{
//  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
//};

struct leptons
{
  TLorentzVector lepton;
//  int pdgId;
//  float dZ;
//  bool passId;
//  bool passVetoId;
//  bool passLooseIso;
//  bool passTightIso;
//  bool passVTightIso;
//  bool passVVTightIso;
};

//struct jets
//{
//  TLorentzVector jet;
//  bool passId;
//  float jetPtJESUp;
//  float jetPtJESDown;
//  float jetEJESUp;
//  float jetEJESDown;
//  float JecUnc;
//};

////lepton highest pt comparator
//struct largest_pt
//{
//  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
//} my_largest_pt;

////jet highest pt comparator
//struct largest_pt_jet
//{
//  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
//} my_largest_pt_jet;


//-------MET

void VLLntupler::Analyze(bool isData, int Option, string outputfilename, string label, string pileupWeightName)
{
 
    cout << "Initializing..." << endl;
    cout << "Label: " <<label<<endl;
    string outfilename = outputfilename;
    if (outfilename == "") outfilename = "VLLNtuple.root";
    TFile *outFile = new TFile(outfilename.c_str(), "RECREATE");    
    TH1F *CutFlow  = new TH1F("CutFlow"   , "CutFlow"   , 5, 0, 5);
    //output TTree
    TTree *outputTree = new TTree("tree", "");
 
    //------------------------
    //declare branch variables
    //------------------------
    bool Trigger;
    float MET;
    float MET_JESUp;
    float MET_JESDn;
    float MET_SF;
    float MET_JESUpSF;
    float MET_JESDnSF;
    float MET_XYCorr;
    float MET_EENoise;
    float MET_HEM;
    float MET_PhiXYCorr;
    float MET_PhiEENoise;
    float MET_PhiHEM;
    float MET_Phi;
    float MET_PhiJESUp;
    float MET_PhiJESDn;
    bool MET_filters;
    bool L1prefiring;
    float MHT;
    //tau
    unsigned int nTaus = 0;//other variable for ID taus
    float tau_M[20];
    float tau_Pt[20];
    float tau_Eta[20];
    float tau_Phi[20];
    float tau_DeltaR[20];
    int   tau_DecayMode[20];
    float tau_Dz[20];
    bool  tau_IsVVVLoose[20];
    bool  tau_IsVVLoose[20];
    bool  tau_IsVLoose[20];
    bool  tau_IsLoose[20];
    bool  tau_IsMedium[20];
    bool  tau_IsTight[20];
    bool  tau_IsVTight[20];
    bool  tau_IsVVTight[20];
    int   tau_GenPartFlav[20];

    //Muon
    unsigned int nMuons = 0;//other variable for ID muons
    float   muon_M[50];
    float   muon_Pt[50];
    float   muon_Eta[50];
    float   muon_Phi[50];
    bool    muon_LooseId[50]; //_looseId

    //Jet
    unsigned int nJets = 0;//other variable for ID muons
    float jet_M[100];
    float jet_Pt[100];
    float jet_PtJESUp[100];
    float jet_PtJESDn[100];
    float jet_Eta[100];
    float jet_Phi[100];
    //float jetBtag[20];
    int   jet_Id[100]; //_jetId: Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
    //int   jetPuId[20];  //_puId : Pilup ID flags with 80X (2016) training 

    //gen tau
    int nGenTau = 0;
    float genTau_Pt[10];
    float genTau_Eta[10];
    float genTau_Phi[10];
    float genTau_M[10];
    
    //GEN Weights
    float gen_Weight;
    int ngenPdfWeights;
    int ngenScaleWeights;
    int ngenPSWeights;
    float gen_PdfWeights[101];
    float gen_ScaleWeights[8];
    float gen_PSWeights[4];

    //load trigger file
    string CMSSWDir  = std::getenv("CMSSW_BASE");
    string trigfile  = CMSSWDir + "/src/TauNanoAODAnalyzer/data/METTriggers_SF.root";
    metTriggerSFFile = TFile::Open(trigfile.c_str());
    if (label=="2016") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Summer16");
    else if (label=="2016UL") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Summer16");
    else if (label=="2016APVUL") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Summer16");
    else if (label=="2017") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall17");
    else if (label=="2017UL") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall17");
    else if (label=="2018") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall18");
    else if (label=="2018UL") metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall18");
    else {metTriggerSFHist = (TH1F*)metTriggerSFFile->Get("trigger_efficiency_Fall18");}    

    //load JEC file
    string JECUncertaintyFile = "";
    if (label == "2016") {
      JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2017") {
      JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2018") {
     JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2016UL") {
     JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Summer19UL16_V7_MC/Summer19UL16_V7_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2016APVUL") {
     JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2017UL") {
     JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Summer19UL17_V5_MC/Summer19UL17_V5_MC_Uncertainty_AK4PFchs.txt";
    } else if (label == "2018UL") {
     JECUncertaintyFile = CMSSWDir + "/src/TauNanoAODAnalyzer/data/JEC/Summer19UL18_V5_MC/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.txt";
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(JECUncertaintyFile.c_str()); 
    //------------------------
    //set branches on big tree
    //------------------------
    //Events info    
    outputTree->Branch("event",             &event,               "event/l");  
    outputTree->Branch("run",               &run,                 "run/i");     
    outputTree->Branch("luminosityBlock",   &luminosityBlock,     "luminosityBlock/i");   

    //GenWeights
    outputTree->Branch("gen_Weight",          &gen_Weight,      "gen_Weight/F");
    outputTree->Branch("ngenPdfWeights",     &ngenPdfWeights, "ngenPdfWeights/i");
    outputTree->Branch( "gen_PdfWeights",      &gen_PdfWeights,  "gen_PdfWeights[ngenPdfWeights]/F");
    outputTree->Branch("ngenScaleWeights",    &ngenScaleWeights, "ngenScaleWeights/i");
    outputTree->Branch( "gen_ScaleWeights",     &gen_ScaleWeights,  "gen_ScaleWeights[ngenScaleWeights]/F");
    outputTree->Branch("ngenPSWeights",        &ngenPSWeights,   "ngenPSWeights/i");
    outputTree->Branch( "gen_PSWeights",         &gen_PSWeights,    "gen_PSWeights[ngenPSWeights]/F"); 

    //Trigger
    outputTree->Branch("Trigger", &Trigger, "Trigger/O"); 

    //EEL1Prefiring
    outputTree->Branch("L1prefiring", &L1prefiring, "L1prefiring/O"); 

    //MET_filters
    outputTree->Branch("MET_filters", &MET_filters, "MET_filters/O"); 

    //MET
    outputTree->Branch("MET", &MET, "MET/F"); 
    outputTree->Branch("MET_JESUp", &MET_JESUp, "MET_JESUp/F"); 
    outputTree->Branch("MET_JESDn", &MET_JESDn, "MET_JESDn/F"); 
    outputTree->Branch("MET_SF", &MET_SF, "MET_SF/F"); 
    outputTree->Branch("MET_JESUpSF", &MET_JESUpSF, "MET_JESUpSF/F"); 
    outputTree->Branch("MET_JESDnSF", &MET_JESDnSF, "MET_JESDnSF/F"); 
    outputTree->Branch("MET_XYCorr" , &MET_XYCorr,  "MET_XYCorr/F"); 
    outputTree->Branch("MET_EENoise", &MET_EENoise, "MET_EENoise/F"); 
    outputTree->Branch("MET_HEM"    , &MET_HEM,     "MET_HEM/F"); 
    outputTree->Branch("MET_Phi", &MET_Phi, "MET_Phi/F"); 
    outputTree->Branch("MET_PhiXYCorr" , &MET_PhiXYCorr,  "MET_PhiXYCorr/F"); 
    outputTree->Branch("MET_PhiEENoise", &MET_PhiEENoise, "MET_PhiEENoise/F"); 
    outputTree->Branch("MET_PhiHEM"    , &MET_PhiHEM,     "MET_PhiHEM/F"); 
    outputTree->Branch("MET_PhiJESUp", &MET_PhiJESUp, "MET_PhiJESUp/F"); 
    outputTree->Branch("MET_PhiJESDn", &MET_PhiJESDn, "MET_PhiJESDn/F"); 
  
    //MHT
    outputTree->Branch("MHT", &MHT, "MHT/F"); 

    //muons
    outputTree->Branch("nMuon",          &nMuons,        "nMuon/i");
    outputTree->Branch("muon_M",          muon_M,          "muon_M[nMuon]/F");
    outputTree->Branch("muon_Pt",         muon_Pt,         "muon_Pt[nMuon]/F");
    outputTree->Branch("muon_Eta",        muon_Eta,        "muon_Eta[nMuon]/F");
    outputTree->Branch("muon_Phi",        muon_Phi,        "muon_Phi[nMuon]/F");
    outputTree->Branch("muon_LooseId",    muon_LooseId,    "muon_LooseId[nMuon]/O");

    //jets
    outputTree->Branch("nJet",  &nJets,                "nJet/i");
    outputTree->Branch("jet_M",        jet_M,       "jet_M[nJet]/F");
    outputTree->Branch("jet_Pt",       jet_Pt,      "jet_Pt[nJet]/F");
    outputTree->Branch("jet_PtJESUp",  jet_PtJESUp, "jet_PtJESUp[nJet]/F");
    outputTree->Branch("jet_PtJESDn",  jet_PtJESDn, "jet_PtJESDn[nJet]/F");
    outputTree->Branch("jet_Eta",      jet_Eta,     "jet_Eta[nJet]/F");
    outputTree->Branch("jet_Phi",      jet_Phi,     "jet_Phi[nJet]/F");
    outputTree->Branch("jet_Id",       jet_Id,      "jet_Id[nJet]/I");

    //tau
    outputTree->Branch("nTau", &nTaus,                  "nTau/i");
    outputTree->Branch("tau_M", tau_M,                    "tau_M[nTau]/F");
    outputTree->Branch("tau_Pt", tau_Pt,                  "tau_Pt[nTau]/F");
    outputTree->Branch("tau_Eta", tau_Eta,                "tau_Eta[nTau]/F");
    outputTree->Branch("tau_Phi", tau_Phi,                "tau_Phi[nTau]/F");
    outputTree->Branch("tau_DeltaR", tau_DeltaR,          "tau_DeltaR[nTau]/F");
    outputTree->Branch("tau_DecayMode", tau_DecayMode,    "tau_DecayMode[nTau]/I");
    outputTree->Branch("tau_IsVVVLoose", tau_IsVVVLoose,"tau_IsVVVLoose[nTau]/O");
    outputTree->Branch("tau_IsVVLoose", tau_IsVVLoose,  "tau_IsVVLoose[nTau]/O");
    outputTree->Branch("tau_IsVLoose", tau_IsVLoose,    "tau_IsVLoose[nTau]/O");
    outputTree->Branch("tau_IsLoose", tau_IsLoose,      "tau_IsLoose[nTau]/O");
    outputTree->Branch("tau_IsMedium", tau_IsMedium,    "tau_IsMedium[nTau]/O");
    outputTree->Branch("tau_IsTight", tau_IsTight,      "tau_IsTight[nTau]/O");
    outputTree->Branch("tau_IsVTight", tau_IsVTight,    "tau_IsVTight[nTau]/O");
    outputTree->Branch("tau_IsVVTight", tau_IsVVTight,  "tau_IsVVTight[nTau]/O");
    outputTree->Branch("tau_Dz",          tau_Dz,            "tau_Dz[nTau]/F");
    outputTree->Branch("tau_GenPartFlav", tau_GenPartFlav,   "tau_GenPartFlav[nTau]/I");
    //gen tau
    outputTree->Branch("nGenTau",          &nGenTau,   "nGenTau/i");
    outputTree->Branch("genTau_Pt",        genTau_Pt,  "genTau_Pt[nGenTau]/F");
    outputTree->Branch("genTau_Eta",       genTau_Eta, "genTau_Eta[nGenTau]/F");
    outputTree->Branch("genTau_Phi",       genTau_Phi, "genTau_Phi[nGenTau]/F");
    outputTree->Branch("genTau_M",         genTau_M,   "genTau_M[nGenTau]/F");
    cout << "Run With Option = " << Option << "\n";
    
    UInt_t NEventsFilled = 0;
   
    //begin loop
    if (fChain == 0) return;
    UInt_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;
    auto start = chrono::steady_clock::now();
    cout << "nentries = " << nentries << "\n";
    //for (UInt_t jentry=0; jentry<1000;jentry++) {
    for (UInt_t jentry=0; jentry<nentries;jentry++) {
      //begin event
      if(jentry % 50000 == 0) {
        cout << "Processing entry " << jentry << endl;
      }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //reset variables
      MET     = -1;
      MET_Phi = -1;
      Trigger=false;
      MET_JESUp = -1;
      MET_JESDn = -1;
      MET_SF = 1;
      MET_JESUpSF = 1;
      MET_JESDnSF = 1;
      MET_XYCorr =-1;
      MET_EENoise=-1;
      MET_HEM    =-1;
      MET_PhiXYCorr=-1;
      MET_PhiEENoise=-1;
      MET_PhiHEM =-1;
      MET_Phi = -1;
      MET_PhiJESUp = -1;
      MET_PhiJESDn = -1;
      MET_filters=false;
      L1prefiring=true;
      nTaus  = 0;
      nMuons  = 0;
      nJets  = 0;
      MHT=0;
      for( int i = 0; i < 20; i++)
     	{
        tau_M[i]  = 0.0;
     	  tau_Pt[i]  = 0.0;
     	  tau_Eta[i] = 0.0;
     	  tau_Phi[i] = 0.0;
        tau_Dz[i]         = 0.0;
     	  tau_DeltaR[i]     = 0.0;
        tau_DecayMode[i]  = 0.0;
        tau_GenPartFlav[i]= 0;
        tau_IsVVVLoose[i]= false;
        tau_IsVVLoose[i] = false;
        tau_IsVLoose[i]  = false;
        tau_IsLoose[i]   = false;
        tau_IsMedium[i]  = false;
        tau_IsTight[i]   = false;
        tau_IsVTight[i]  = false;
        tau_IsVVTight[i] = false;
     	}

      for( int i = 0; i < 50; i++)
      {
        muon_M[i]  = 0.0;
        muon_Pt[i]  = 0.0;
        muon_Eta[i] = 0.0;
        muon_Phi[i] = 0.0;
        muon_LooseId[i]     = false;
      }

      for( int i = 0; i < 100; i++)
      {
        jet_M[i]   = 0.0;
        jet_Pt[i]  = 0.0;
        jet_Eta[i] = 0.0;
        jet_Phi[i] = 0.0;
        jet_Id[i]  = 0.0;
      }

      nGenTau = 0;
      for( int i = 0; i < 5; i++)
      {
	      genTau_Pt[i]  = 0.0;
	      genTau_Eta[i] = 0.0;
	      genTau_Phi[i] = 0.0;
        genTau_M[i] = 0.0;
	    }

      //GenWeights (if MC)
      if (!isData)
      {
         gen_Weight=genWeight;
         ngenPdfWeights=nLHEPdfWeight;
         ngenScaleWeights=nLHEScaleWeight;
         ngenPSWeights=nPSWeight;
         for( int i = 0; i < ngenPdfWeights; i++)
         {
           gen_PdfWeights[i]=LHEPdfWeight[i];
         }
         for( int i = 0; i < ngenScaleWeights; i++)
         {
           gen_ScaleWeights[i]=LHEScaleWeight[i];
         }
         for( int i = 0; i < ngenPSWeights; i++)
         {
           gen_PSWeights[i]=PSWeight[i];
         }  
      }      
      //*************************
      //Trigger
      //*************************
      if ( (label=="2016"  ||  label=="2016UL" || label=="2016APVUL") && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==true) )  Trigger=true; 
      if ( (label=="2017" || label=="2018" || label=="2017UL" || label=="2018UL") && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight==true || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight==true || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60==true) )  Trigger=true;
      CutFlow->Fill(0);
      if (isData && Trigger==false) continue; //Trigger cut on the data only

      //*************************
      //MET filters
      //*************************
      bool passfilters = false;
      if (Flag_goodVertices==true && Flag_globalSuperTightHalo2016Filter==true && Flag_HBHENoiseFilter==true && Flag_EcalDeadCellTriggerPrimitiveFilter==true && Flag_BadPFMuonFilter==true && Flag_eeBadScFilter==true) passfilters=true;
      if (isData)
      {
        if(label=="2016UL" && passfilters == true && Flag_BadPFMuonDzFilter == true)  MET_filters=true;
        else if (label=="2016APVUL" && passfilters == true && Flag_BadPFMuonDzFilter == true)  MET_filters=true;
        else if (label=="2017UL" && passfilters == true && Flag_ecalBadCalibFilter==true && Flag_BadPFMuonDzFilter == true) MET_filters=true;
        else if (label=="2018UL" && passfilters == true && Flag_ecalBadCalibFilter==true && Flag_BadPFMuonDzFilter == true) MET_filters=true;
        else if (label=="2016" && passfilters == true)  MET_filters=true;
        else if (label=="2017" && passfilters == true && Flag_ecalBadCalibFilter==true) MET_filters=true;
        else if (label=="2018" && passfilters == true && Flag_ecalBadCalibFilter==true) MET_filters=true;
      }
      else{
        if(label=="2016" && passfilters == true)  MET_filters=true;
        else if (label=="2017" && passfilters == true && Flag_ecalBadCalibFilter==true) MET_filters=true;
        else if (label=="2018" && passfilters == true && Flag_ecalBadCalibFilter==true) MET_filters=true;         
        else if (label=="2016UL" && passfilters == true && Flag_BadPFMuonDzFilter == true)  MET_filters=true;
        else if (label=="2016APVUL" && passfilters == true && Flag_BadPFMuonDzFilter == true)  MET_filters=true;
        else if (label=="2017UL" && passfilters == true && Flag_ecalBadCalibFilter==true && Flag_BadPFMuonDzFilter == true) MET_filters=true;
        else if (label=="2018UL" && passfilters == true && Flag_ecalBadCalibFilter==true && Flag_BadPFMuonDzFilter == true) MET_filters=true;   
      }

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************
      std::vector<leptons> Leptons;      
      //********************************
      //Leptons (for overlaps)
      //********************************
      for(int i = 0; i < nMuon;  i++){
        if(Muon_pt[i] < 25) continue;
        if(abs(Muon_eta[i]) > 2.4) continue;
        if(!Muon_looseId[i]) continue;
        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (EventAnalyzer::deltaR(Muon_eta[i],Muon_phi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i], Muon_phi[i],MU_MASS);
//        tmpMuon.pdgId = 13 * -1 * Muon_charge[i];
//        tmpMuon.dZ = Muon_dz[i];
//        tmpMuon.passId = Muon_tightId[i];
//        if ( Muon_pfIsoId[i] >= 2 ) tmpMuon.passLooseIso  = true;
//        if ( Muon_pfIsoId[i] >= 4 ) tmpMuon.passTightIso  = true;
//        if ( Muon_pfIsoId[i] >= 5 ) tmpMuon.passVTightIso = true;
//        if ( Muon_pfIsoId[i] >= 6 ) tmpMuon.passVVTightIso= true;
//        tmpMuon.passVetoId = false;
        Leptons.push_back(tmpMuon);       
      }
      for( int i = 0; i < nElectron; i++ )
      {
          if (!Electron_mvaFall17V2Iso_WPL[i]) continue; //LooseID
          if(Electron_pt[i] < 35) continue;
          if(fabs(Electron_eta[i]) > 2.4) continue;
          //remove overlaps
          bool overlap = false;
          for(auto& lep : Leptons)
          {
            if (EventAnalyzer::deltaR(Electron_eta[i],Electron_phi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
          }
          if(overlap) continue;
          leptons tmpElectron;
          tmpElectron.lepton.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i], Electron_phi[i], ELE_MASS);
//          tmpElectron.pdgId = 11 * -1 * Electron_charge[i];
//          tmpElectron.dZ = Electron_dz[i];
//          tmpElectron.passId     = Electron_mvaFall17V2Iso_WP80[i]; //TightID
//          tmpElectron.passVetoId = false;
          Leptons.push_back(tmpElectron);
      }


      //*************************
      //MET variables
      //*************************
      MET       = MET_pt;
      MET_Phi   = MET_phi;
      MET_JESUp = MET_pt;
      MET_JESDn = MET_pt;
      if(!isData) MET_SF    = EventAnalyzer::GetMETTriggerSF(MET_pt,metTriggerSFHist);
      std::pair<double,double> corrected_met;
      if (label=="2016") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016", !isData, PV_npvs,false);
      else if (label=="2017") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2017", !isData, PV_npvs,false);
      else if (label=="2018") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2018", !isData, PV_npvs,false);
      else if (label=="2016UL") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016nonAPV", !isData, PV_npvs, true);
      else if (label=="2016APVUL") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2016APV", !isData, PV_npvs, true);
      else if (label=="2017UL") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2017", !isData, PV_npvs, true);
      else if (label=="2018UL") corrected_met = EventAnalyzer::METXYCorr_Met_MetPhi(MET_pt, MET_phi, run, "2018", !isData, PV_npvs, true);
      MET_XYCorr     = corrected_met.first;
      MET_EENoise    = MET_XYCorr;
      MET_HEM        = MET_XYCorr;
      MET_PhiXYCorr  = corrected_met.second;
      MET_PhiEENoise = MET_PhiXYCorr;
      MET_PhiHEM     = MET_PhiXYCorr;

      //-----------------------------------------------
      //MET/Jets
      //-----------------------------------------------
      float MetXCorr_JESUp = 0.;
      float MetYCorr_JESUp = 0.;
      float MetXCorr_JESDown = 0.;
      float MetYCorr_JESDown = 0.;
      float MetXCorr_HEM = 0.;
      float MetYCorr_HEM = 0.;
      float MetXCorr_EENoise = 0.;
      float MetYCorr_EENoise = 0.;
      float MHT_X=0.;
      float MHT_Y=0.;
      for(int i = 0; i < nJet; i++)
      {
          //For MHT calculation
          TLorentzVector theJet = EventAnalyzer::makeTLorentzVectorPtEtaPhiM( Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i] );


          //offline jetid
          bool offlinejetid=false;
          if(Jet_chEmEF[i] < 0.99 && Jet_chHEF[i] > 0 && Jet_neEmEF[i] < 0.99 && Jet_neHEF[i] < 0.90 && Jet_nConstituents[i] > 1) offlinejetid=true;
          //compute MHT ALA NANOAOD PRODUCER (noID)
          if(fabs(theJet.Eta()) < 5.0 && theJet.Pt() > 20 && offlinejetid==true) //&& Jet_jetId[i]>=2
          {
            MHT_X += -theJet.Px();
            MHT_Y += -theJet.Py();
          }
          //------------------------------------------------------------
          //exclude selected muons and electrons from the jet collection
          //------------------------------------------------------------
          double deltaR = -1;
          for(auto& lep : Leptons){
            double thisDR = EventAnalyzer::deltaR(Jet_eta[i],Jet_phi[i],lep.lepton.Eta(),lep.lepton.Phi());
            if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
          }
          if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton
          double jetCorrPt = Jet_pt[i];
          double jetCorrM  = Jet_mass[i];
          TLorentzVector thisJet = EventAnalyzer::makeTLorentzVectorPtEtaPhiM( jetCorrPt, Jet_eta[i], Jet_phi[i], jetCorrM );
          //Era dependent special corrections to MET (HEM=Remove jets in affected region, EENoise:Remove noise jets, EEprefiring: Tag events with jets in prefiring region)
          if (thisJet.Eta()>-3.0 && thisJet.Eta()<-1.3 && thisJet.Phi() >-1.57 && thisJet.Phi() <-0.87 && (label == "2018" || label == "2018UL") )
          {
            MetXCorr_HEM += thisJet.Px();
            MetYCorr_HEM += thisJet.Py();
          }
          if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && (label == "2017" || label == "2017UL" )  )
          {
            MetXCorr_EENoise += thisJet.Px();
            MetYCorr_EENoise += thisJet.Py();
          }
          if (fabs(thisJet.Eta())> 2.25 && fabs(thisJet.Eta())<3.0 && thisJet.Pt() > 100 )
          {
            L1prefiring = false;
          }
          if (fabs(thisJet.Eta()) >= 3.0) continue;
          //JES variations
          if (!isData)
          {
              jecUnc->setJetEta(Jet_eta[i]);
              jecUnc->setJetPt(Jet_pt[i]);
              double unc = jecUnc->getUncertainty(true);
              TLorentzVector thisJetJESUp   = EventAnalyzer::makeTLorentzVector(Jet_pt[i]*(1+unc), Jet_eta[i], Jet_phi[i], thisJet.E()*(1+unc) );
              TLorentzVector thisJetJESDown = EventAnalyzer::makeTLorentzVector(Jet_pt[i]*(1-unc), Jet_eta[i], Jet_phi[i], thisJet.E()*(1-unc) );
              if (thisJetJESUp.Pt() > 10)//15 && Jet_neEmEF[i]<0.9)
              {
                MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
                MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
              }
              if (thisJetJESDown.Pt() > 10)//15 && Jet_neEmEF[i]<0.9)
              {
                MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
                MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
              }
          }
      }

      //Fill MHT
      MHT  = sqrt( pow(MHT_X,2) + pow(MHT_Y,2) );

      //Fill era dependent special corrections to MET + JES variation effects to MET: Nominal branch is called EENoise (Christina's)
      TLorentzVector PFMET = EventAnalyzer::makeTLorentzVectorPtEtaPhiM(MET_XYCorr, 0, MET_PhiXYCorr, 0);
      if(!isData)
      {
        //JES up
        float PFMetXJESUp          = PFMET.Px() + MetXCorr_JESUp;
        float PFMetYJESUp          = PFMET.Py() + MetYCorr_JESUp;
        MET_JESUp                  = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
        MET_PhiJESUp               = atan(PFMetYJESUp/PFMetXJESUp);
        if  (PFMetXJESUp < 0.0) MET_PhiJESUp = EventAnalyzer::deltaPhi(TMath::Pi() + MET_PhiJESUp,0.0);
        MET_JESUpSF                = EventAnalyzer::GetMETTriggerSF(MET_JESUp,metTriggerSFHist)/MET_SF;
        //JES down
        float PFMetXJESDown        = PFMET.Px() + MetXCorr_JESDown;
        float PFMetYJESDown        = PFMET.Py() + MetYCorr_JESDown;
        MET_JESDn                  =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
        MET_PhiJESDn               = atan(PFMetYJESDown/PFMetXJESDown);
        if  (PFMetXJESUp < 0.0) MET_PhiJESDn = EventAnalyzer::deltaPhi(TMath::Pi() + MET_PhiJESDn,0.0);
        MET_JESDnSF                = EventAnalyzer::GetMETTriggerSF(MET_JESDn,metTriggerSFHist)/MET_SF;
      }
      //HEM
      float PFMetXHEM            = PFMET.Px() + MetXCorr_HEM;
      float PFMetYHEM            = PFMET.Py() + MetYCorr_HEM;
      MET_HEM         = sqrt( pow(PFMetXHEM,2) + pow(PFMetYHEM,2) );
      MET_PhiHEM      = atan(PFMetYHEM/PFMetXHEM);
      if  (PFMetXHEM < 0.0) MET_PhiHEM = EventAnalyzer::deltaPhi(TMath::Pi() + MET_PhiHEM,0.0);
      //EENoise
      float PFMetXEENoise        = PFMET.Px() + MetXCorr_EENoise;
      float PFMetYEENoise        = PFMET.Py() + MetYCorr_EENoise;
      MET_EENoise     = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
      MET_PhiEENoise  = atan(PFMetYEENoise/PFMetXEENoise);
      if  (PFMetXEENoise < 0.0) MET_PhiEENoise = EventAnalyzer::deltaPhi(TMath::Pi() + MET_PhiEENoise,0.0);
      if (isData && MET_EENoise < 200) continue; //Offline MET cut only on data
      CutFlow->Fill(1);

      //********************************
      //Tau
      //********************************
      for(int i = 0; i < nTau;  i++){
        if(Tau_pt[i] < 20) continue;
        if(abs(Tau_eta[i]) > 2.3) continue;
        if(abs(Tau_dz[i]) > 0.2) continue;
        if(Tau_decayMode[i]==5 || Tau_decayMode[i]==6) continue; 
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 1  && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsVVVLoose[nTaus] = true; 
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 2  && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsVVLoose[nTaus]  = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 4  && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsVLoose[nTaus]   = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 8  && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsLoose[nTaus]    = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 16 && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsMedium[nTaus]   = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 32 && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsTight[nTaus]    = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >= 64 && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsVTight[nTaus]   = true;
        if ( Tau_idDeepTau2017v2p1VSjet[i] >=128 && Tau_idDeepTau2017v2p1VSe[i] >= 2 && Tau_idDeepTau2017v2p1VSmu[i] >= 8 ) tau_IsVVTight[nTaus]  = true;
        tau_M[nTaus]            = Tau_mass[i];
        tau_Pt[nTaus]           = Tau_pt[i];
        tau_Eta[nTaus]          = Tau_eta[i];
        tau_Phi[nTaus]          = Tau_phi[i];
        tau_Dz[nTaus]           = Tau_dz[i];
        tau_GenPartFlav[nTaus]  = Tau_genPartFlav[i];
        tau_DecayMode[nTaus]    = Tau_decayMode[i];
        nTaus++;
        TLorentzVector thisTau; thisTau.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
      }
      if (!isData) 
      {
            //********************************
            //Gen Taus
            //********************************
            for( int i = 0; i < nGenPart; i++)
            {
              if( abs(GenPart_pdgId[i]) ==  15 && GenPart_status[i] == 2  && abs(GenPart_pdgId[FindMotherIndex(i)]) == 5000002) //24 //5000002
                {
                  genTau_Pt[nGenTau]  = GenPart_pt[i];
                  genTau_Eta[nGenTau] = GenPart_eta[i];
                  genTau_Phi[nGenTau] = GenPart_phi[i];
                  genTau_M[nGenTau]   = 1.77686;
                  nGenTau++;
                }
            }
            //Reco Taus to Gen taus matching
            for (int i = 0; i < nTaus; i++)
            {
              float minDeltaR = 999.;
              TLorentzVector thisTau; thisTau.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
              for( int j = 0; j < nGenTau; j++ )
                {
                  TLorentzVector thisGenTau; thisGenTau.SetPtEtaPhiM(genTau_Pt[j], genTau_Eta[j], genTau_Phi[j], genTau_M[nGenTau]);
                  if ( thisTau.DeltaR( thisGenTau ) < minDeltaR )
                  {
                    minDeltaR = thisTau.DeltaR( thisGenTau );
                  }
                }
              tau_DeltaR[i] = minDeltaR;
            }
      }   
      //*************************
      //Muons for ClusterVeto
      //*************************
      for(int i = 0; i < nMuon;  i++){
        if(abs(Muon_eta[i]) > 3.0) continue;
        muon_M[nMuons]        = MU_MASS;
        muon_Pt[nMuons]       = Muon_pt[i];
        muon_Eta[nMuons]      = Muon_eta[i];
        muon_Phi[nMuons]      = Muon_phi[i];
        muon_LooseId[nMuons]  = Muon_looseId[i];
        nMuons++;
      }
      //*************************
      //Jets for ClusterVeto
      //*************************
      for(int i = 0; i < nJet;  i++){
        if(abs(Jet_eta[i]) > 3.0) continue;
        jet_M[nJets]          = Jet_mass[i];
        jet_Pt[nJets]         = Jet_pt[i];
        jet_PtJESUp[nJets]    = Jet_pt[i];
        jet_PtJESDn[nJets]    = Jet_pt[i];
        if(!isData)
        {
           jecUnc->setJetEta(Jet_eta[i]);
           jecUnc->setJetPt(Jet_pt[i]);
           double unc = jecUnc->getUncertainty(true);
           jet_PtJESUp[nJets]    = Jet_pt[i]*(1+unc);
           jet_PtJESDn[nJets]    = Jet_pt[i]*(1-unc);
        }
        jet_Eta[nJets]        = Jet_eta[i];
        jet_Phi[nJets]        = Jet_phi[i];
        jet_Id[nJets]         = Jet_jetId[i];
        nJets++;
      }
      //-----------------------------------------------
      //Event selection and cutflow histogram 
      //-----------------------------------------------
      //Bin 0: Total number of events
      //CutFlow->Fill(0);
      //Bin 1: MET Trigger + MET cut
      //if ( Trigger==false ) continue;  
      //if ( MET < 200 ) continue;
      //CutFlow->Fill(1);
      //***********************
      //Fill Event
      //***********************
      NEventsFilled++;
      outputTree->Fill();      
    }//end of event loop
    auto end = chrono::steady_clock::now();
    auto diff = end - start;
    cout << "Filled Total of " << NEventsFilled << " Events\n";
    cout << "Writing output trees..." << endl;
    cout <<"Total time: " << chrono::duration <double, milli> (diff).count() << " ms" << endl;
    //outFile->cd();
    outFile->Write();
    outFile->Close();

}



