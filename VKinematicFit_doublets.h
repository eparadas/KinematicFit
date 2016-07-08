#ifndef VKinematicFit_doublets_h
#define VKinematicFit_doublets_h


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMatrixD.h>
// #include "VDecompLU.h"
#include <TDecompLU.h>
#include <TApplication.h>

using namespace std;

class VKinematicFit_doublets
{
private:
  float initialPt[8], correctedPt[8], initialEta[8], initialPhi[8], Corrections[11], zeroCorrections[8];
  float s_lamda[3];
  int nIterations;
  int countedIterations;
  double Dif_Chi2;
  double final_chi2;
  float initialDoublet_1Mass[2], initialDoublet_2Mass[2], correctedDoublet_1Mass[2], correctedDoublet_2Mass[2], correctedQuartetmass[2];
  
public:
  VKinematicFit_doublets(float *f_pt, float *f_eta, float *f_phi);
  void initialize();
  void setIterations(int iters);
  void setChi2Dif(double dif);
  void getCorrectedPt(float *new_pt);
  void getMass_from_Selection(float *masses_1, float *masses_2);
  void getMass_after_Correction(float *masses_1, float *masses_2);
  void getMass_for_Quartets(float *masses);
  float getChi2Value();
  int ApplyFit();
};

#endif