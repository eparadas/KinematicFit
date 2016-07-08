//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov  1 11:59:03 2013 by ROOT version 5.34/05
// from TTree events/events
// found on file: /mnt/storage/8TeV/FastSimulation/CMG_1200_300_120/flatTree_1200_300_120.root
//////////////////////////////////////////////////////////

#ifndef flat_Class_h
#define flat_Class_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxrun = 1;
const Int_t kMaxevt = 1;
const Int_t kMaxlumi = 1;
const Int_t kMaxnVtx = 1;
const Int_t kMaxpvx = 1;
const Int_t kMaxpvy = 1;
const Int_t kMaxpvz = 1;
const Int_t kMaxindex2J = 1;
const Int_t kMaxindex4J = 1;
const Int_t kMaxrho = 1;
const Int_t kMaxmetSig = 1;
const Int_t kMaxdPhi4j = 1;
const Int_t kMaxht = 1;
const Int_t kMaxm8j = 1;
const Int_t kMaxm4jAve = 1;
const Int_t kMaxm2jAve = 1;
const Int_t kMaxm2jEst = 1;
const Int_t kMaxm2jSigma = 1;
const Int_t kMaxcosThetaStar = 1;
const Int_t kMaxm4jBalance = 1;
const Int_t kMaxm2j = 1;
const Int_t kMaxm4j = 1;
const Int_t kMaxht4j = 1;
const Int_t kMaxeta4j = 1;
const Int_t kMaxpt4j = 1;
const Int_t kMaxdR2jAll = 1;
const Int_t kMaxsimPU = 1;
const Int_t kMaxm4jAveGEN = 1;
const Int_t kMaxm2jAveGEN = 1;
const Int_t kMaxm2jEstGEN = 1;
const Int_t kMaxm4jAveParton = 1;
const Int_t kMaxm2jAveParton = 1;
const Int_t kMaxm2jEstParton = 1;

class flat_Class {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumi;
   Int_t           nvtx;
   Float_t         pvx;
   Float_t         pvy;
   Float_t         pvz;
   Int_t           index2J[4][2];
   Int_t           index4J[2][4];
   Float_t         rho;
   Float_t         metSig;
   Float_t         dphi4j;
   Float_t         ht;
   Float_t         m8j;
   Float_t         m4jAve;
   Float_t         m2jAve;
   Float_t         m2jEst;
   Float_t         m2jSig;
   Float_t         cosThetaStar;
   Float_t         m4jBalance;
   Float_t         m2j[4];
   Float_t         m4j[2];
   Float_t         ht4j[2];
   Float_t         eta4j[2];
   Float_t         pt4j[2];
   Float_t         dR2jAll[28];
   Float_t         total_ht;
   vector<float>   *pt;
   vector<float>   *jec;
   vector<float>   *unc;
   vector<float>   *beta;
   vector<float>   *eta;
   vector<float>   *phi;
   vector<float>   *mass;
   vector<float>   *chf;
   vector<float>   *nhf;
   vector<float>   *phf;
   vector<float>   *muf;
   vector<float>   *elf;
   vector<int>     *bTagIdx;
   vector<int>     *etaIdx;
   vector<float>   *btag;
   vector<float>   *puMva;
   vector<int>     *puId;
   vector<float>   *ptD;
   vector<float>   *ptD_QC;
   vector<float>   *AxisMinor;
   vector<float>   *AxisMajor;
   vector<float>   *AxisMinor_QC;
   vector<float>   *AxisMajor_QC;
   vector<float>   *pull;
   vector<float>   *pull_QC;
   vector<float>   *jetR;
   vector<float>   *jetRChg_QC;
   Int_t           pu;
   Float_t         m4jAveGEN;
   Float_t         m2jAveGEN;
   Float_t         m2jEstGEN;
   Float_t         m4jAveParton;
   Float_t         m2jAveParton;
   Float_t         m2jEstParton;
   vector<int>     *partonId;
   vector<int>     *partonSt;
   vector<float>   *partonPt;
   vector<float>   *partonEta;
   vector<float>   *partonPhi;
   vector<float>   *partonE;
   vector<float>   *genjetPt;
   vector<float>   *genjetEta;
   vector<float>   *genjetPhi;
   vector<float>   *genjetE;

   // List of branches
   TBranch        *b_run_;   //!
   TBranch        *b_evt_;   //!
   TBranch        *b_lumi_;   //!
   TBranch        *b_nVtx_;   //!
   TBranch        *b_pvx_;   //!
   TBranch        *b_pvy_;   //!
   TBranch        *b_pvz_;   //!
   TBranch        *b_index2J_;   //!
   TBranch        *b_index4J_;   //!
   TBranch        *b_rho_;   //!
   TBranch        *b_metSig_;   //!
   TBranch        *b_dPhi4j_;   //!
   TBranch        *b_ht_;   //!
   TBranch        *b_m8j_;   //!
   TBranch        *b_m4jAve_;   //!
   TBranch        *b_m2jAve_;   //!
   TBranch        *b_m2jEst_;   //!
   TBranch        *b_m2jSigma_;   //!
   TBranch        *b_cosThetaStar_;   //!
   TBranch        *b_m4jBalance_;   //!
   TBranch        *b_m2j_;   //!
   TBranch        *b_m4j_;   //!
   TBranch        *b_ht4j_;   //!
   TBranch        *b_eta4j_;   //!
   TBranch        *b_pt4j_;   //!
   TBranch        *b_dR2jAll_;   //!
   TBranch        *b_total_ht;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_jec;   //!
   TBranch        *b_unc;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_chf;   //!
   TBranch        *b_nhf;   //!
   TBranch        *b_phf;   //!
   TBranch        *b_muf;   //!
   TBranch        *b_elf;   //!
   TBranch        *b_bTagIdx;   //!
   TBranch        *b_etaIdx;   //!
   TBranch        *b_btag;   //!
   TBranch        *b_puMva;   //!
   TBranch        *b_puId;   //!
   TBranch        *b_ptD;   //!
   TBranch        *b_ptD_QC;   //!
   TBranch        *b_AxisMinor;   //!
   TBranch        *b_AxisMajor;   //!
   TBranch        *b_AxisMinor_QC;   //!
   TBranch        *b_AxisMajor_QC;   //!
   TBranch        *b_pull;   //!
   TBranch        *b_pull_QC;   //!
   TBranch        *b_jetR;   //!
   TBranch        *b_jetRChg_QC;   //!
   TBranch        *b_simPU_;   //!
   TBranch        *b_m4jAveGEN_;   //!
   TBranch        *b_m2jAveGEN_;   //!
   TBranch        *b_m2jEstGEN_;   //!
   TBranch        *b_m4jAveParton_;   //!
   TBranch        *b_m2jAveParton_;   //!
   TBranch        *b_m2jEstParton_;   //!
   TBranch        *b_partonId;   //!
   TBranch        *b_partonSt;   //!
   TBranch        *b_partonPt;   //!
   TBranch        *b_partonEta;   //!
   TBranch        *b_partonPhi;   //!
   TBranch        *b_partonE;   //!
   TBranch        *b_genjetPt;   //!
   TBranch        *b_genjetEta;   //!
   TBranch        *b_genjetPhi;   //!
   TBranch        *b_genjetE;   //!

   flat_Class(TTree *tree=0);
   virtual ~flat_Class();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef flat_Class_cxx
flat_Class::flat_Class(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/storage/8TeV/FastSimulation/CMG_1200_300_120/flatTree_1200_300_120.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mnt/storage/8TeV/FastSimulation/CMG_1200_300_120/flatTree_1200_300_120.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/mnt/storage/8TeV/FastSimulation/CMG_1200_300_120/flatTree_1200_300_120.root:/Multijets");
      dir->GetObject("events",tree);

   }
   Init(tree);
}

flat_Class::~flat_Class()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t flat_Class::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t flat_Class::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void flat_Class::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pt = 0;
   jec = 0;
   unc = 0;
   beta = 0;
   eta = 0;
   phi = 0;
   mass = 0;
   chf = 0;
   nhf = 0;
   phf = 0;
   muf = 0;
   elf = 0;
   bTagIdx = 0;
   etaIdx = 0;
   btag = 0;
   puMva = 0;
   puId = 0;
   ptD = 0;
   ptD_QC = 0;
   AxisMinor = 0;
   AxisMajor = 0;
   AxisMinor_QC = 0;
   AxisMajor_QC = 0;
   pull = 0;
   pull_QC = 0;
   jetR = 0;
   jetRChg_QC = 0;
   partonId = 0;
   partonSt = 0;
   partonPt = 0;
   partonEta = 0;
   partonPhi = 0;
   partonE = 0;
   genjetPt = 0;
   genjetEta = 0;
   genjetPhi = 0;
   genjetE = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNo", &runNo, &b_run_);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evt_);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nVtx_);
   fChain->SetBranchAddress("pvx", &pvx, &b_pvx_);
   fChain->SetBranchAddress("pvy", &pvy, &b_pvy_);
   fChain->SetBranchAddress("pvz", &pvz, &b_pvz_);
   fChain->SetBranchAddress("index2J", index2J, &b_index2J_);
   fChain->SetBranchAddress("index4J", index4J, &b_index4J_);
   fChain->SetBranchAddress("rho", &rho, &b_rho_);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig_);
   fChain->SetBranchAddress("dphi4j", &dphi4j, &b_dPhi4j_);
   fChain->SetBranchAddress("ht", &ht, &b_ht_);
   fChain->SetBranchAddress("m8j", &m8j, &b_m8j_);
   fChain->SetBranchAddress("m4jAve", &m4jAve, &b_m4jAve_);
   fChain->SetBranchAddress("m2jAve", &m2jAve, &b_m2jAve_);
   fChain->SetBranchAddress("m2jEst", &m2jEst, &b_m2jEst_);
   fChain->SetBranchAddress("m2jSig", &m2jSig, &b_m2jSigma_);
   fChain->SetBranchAddress("cosThetaStar", &cosThetaStar, &b_cosThetaStar_);
   fChain->SetBranchAddress("m4jBalance", &m4jBalance, &b_m4jBalance_);
   fChain->SetBranchAddress("m2j", m2j, &b_m2j_);
   fChain->SetBranchAddress("m4j", m4j, &b_m4j_);
   fChain->SetBranchAddress("ht4j", ht4j, &b_ht4j_);
   fChain->SetBranchAddress("eta4j", eta4j, &b_eta4j_);
   fChain->SetBranchAddress("pt4j", pt4j, &b_pt4j_);
   fChain->SetBranchAddress("dR2jAll", dR2jAll, &b_dR2jAll_);
   fChain->SetBranchAddress("total_ht", &total_ht, &b_total_ht);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("jec", &jec, &b_jec);
   fChain->SetBranchAddress("unc", &unc, &b_unc);
   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("mass", &mass, &b_mass);
   fChain->SetBranchAddress("chf", &chf, &b_chf);
   fChain->SetBranchAddress("nhf", &nhf, &b_nhf);
   fChain->SetBranchAddress("phf", &phf, &b_phf);
   fChain->SetBranchAddress("muf", &muf, &b_muf);
   fChain->SetBranchAddress("elf", &elf, &b_elf);
   fChain->SetBranchAddress("bTagIdx", &bTagIdx, &b_bTagIdx);
   fChain->SetBranchAddress("etaIdx", &etaIdx, &b_etaIdx);
   fChain->SetBranchAddress("btag", &btag, &b_btag);
   fChain->SetBranchAddress("puMva", &puMva, &b_puMva);
   fChain->SetBranchAddress("puId", &puId, &b_puId);
   fChain->SetBranchAddress("ptD", &ptD, &b_ptD);
   fChain->SetBranchAddress("ptD_QC", &ptD_QC, &b_ptD_QC);
   fChain->SetBranchAddress("AxisMinor", &AxisMinor, &b_AxisMinor);
   fChain->SetBranchAddress("AxisMajor", &AxisMajor, &b_AxisMajor);
   fChain->SetBranchAddress("AxisMinor_QC", &AxisMinor_QC, &b_AxisMinor_QC);
   fChain->SetBranchAddress("AxisMajor_QC", &AxisMajor_QC, &b_AxisMajor_QC);
   fChain->SetBranchAddress("pull", &pull, &b_pull);
   fChain->SetBranchAddress("pull_QC", &pull_QC, &b_pull_QC);
   fChain->SetBranchAddress("jetR", &jetR, &b_jetR);
   fChain->SetBranchAddress("jetRChg_QC", &jetRChg_QC, &b_jetRChg_QC);
   fChain->SetBranchAddress("pu", &pu, &b_simPU_);
   fChain->SetBranchAddress("m4jAveGEN", &m4jAveGEN, &b_m4jAveGEN_);
   fChain->SetBranchAddress("m2jAveGEN", &m2jAveGEN, &b_m2jAveGEN_);
   fChain->SetBranchAddress("m2jEstGEN", &m2jEstGEN, &b_m2jEstGEN_);
   fChain->SetBranchAddress("m4jAveParton", &m4jAveParton, &b_m4jAveParton_);
   fChain->SetBranchAddress("m2jAveParton", &m2jAveParton, &b_m2jAveParton_);
   fChain->SetBranchAddress("m2jEstParton", &m2jEstParton, &b_m2jEstParton_);
   fChain->SetBranchAddress("partonId", &partonId, &b_partonId);
   fChain->SetBranchAddress("partonSt", &partonSt, &b_partonSt);
   fChain->SetBranchAddress("partonPt", &partonPt, &b_partonPt);
   fChain->SetBranchAddress("partonEta", &partonEta, &b_partonEta);
   fChain->SetBranchAddress("partonPhi", &partonPhi, &b_partonPhi);
   fChain->SetBranchAddress("partonE", &partonE, &b_partonE);
   fChain->SetBranchAddress("genjetPt", &genjetPt, &b_genjetPt);
   fChain->SetBranchAddress("genjetEta", &genjetEta, &b_genjetEta);
   fChain->SetBranchAddress("genjetPhi", &genjetPhi, &b_genjetPhi);
   fChain->SetBranchAddress("genjetE", &genjetE, &b_genjetE);
   Notify();
}

Bool_t flat_Class::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void flat_Class::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t flat_Class::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef flat_Class_cxx
