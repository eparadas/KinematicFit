#ifndef VKinematicFit_quartets_h
#define VKinematicFit_quartets_h


#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMatrixD.h>
#include <TDecompLU.h>
// #include "VDecompLU.h"
#include <TApplication.h>

using namespace std;

class VKinematicFit_quartets
{
private:
  float initialPt[8], correctedPt[8], initialEta[8], initialPhi[8], Corrections[9], zeroCorrections[8];
  float s_lamda;
  int nIterations;
  int countedIterations;
  double Dif_Chi2;
  double final_chi2;
  float initialColoronMass[2], correctedColoronMass[2];
  
public:
  VKinematicFit_quartets(float *f_pt, float *f_eta, float *f_phi);
  void initialize();
  void setIterations(int iters);
  void setChi2Dif(double dif);
  void getCorrectedPt(float *new_pt);
  void getMass_from_Selection(float *masses);
  void getMass_after_Correction(float *masses);
  float getChi2Value();
  int ApplyFit();
};

#endif