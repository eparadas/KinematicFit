#include "./VKinematicFit_quartets.h"
#include "./nVDecompLU.cpp"
#include "./KinFit_class_with_constraints_quadrupal_only_quartets.cpp"

VKinematicFit_quartets::VKinematicFit_quartets(float *f_pt, float *f_eta, float *f_phi)
{
  for(int i=0; i<8; i++)
  {
    initialPt[i] = f_pt[i];
    initialEta[i]= f_eta[i];
    initialPhi[i]= f_phi[i];
    Corrections[i] = 0;
    zeroCorrections[i] = 0;
  }
  Corrections[8] = 0;
}

void VKinematicFit_quartets::initialize()
{
  s_lamda = 1.0;
  nIterations = 200;
  countedIterations = 0;
  Dif_Chi2 = 1e-12;
  final_chi2 = 100000000;
}

int VKinematicFit_quartets::ApplyFit()
{
  initialize();
  bool decomp_flag = false;
  Pts Combination_unchanged( initialPt, initialEta, initialPhi, s_lamda );
			
  int iter = 0;
  double prev_chi, dif_chi;
  do{
    iter++;
    countedIterations++;
    prev_chi = final_chi2;
    for (int k=0; k<8; k++)	
      initialPt[k] += Corrections[k];
      
    s_lamda += Corrections[8];

    
    Pts Combination( initialPt, initialEta, initialPhi, s_lamda );
    
    double const_matrix[9];
    for (int l=0; l<9; l++)
      const_matrix[l] = (Combination.Derived_pt(l,9) + Combination.Derived_quartet_term(l,9)+ Combination.Derived_q_lamda_term(l,9));
   
    TMatrixD final_matrix(9,9);
    for (int l=0; l<9; l++)
    {
      for (int m=0; m<9; m++)
	final_matrix(l,m) = (Combination.Derived_pt(l,m) + Combination.Derived_quartet_term(l,m)+ Combination.Derived_q_lamda_term(l,m));
    }

// 		if ( (! nVDecompLU(final_matrix)) || (final_matrix.Determinant() == 0) )
// 		{
// 			decomp_flag = true;
// 			break;
//     }
		TDecompLU lu(final_matrix);
		try
		{
			lu.Invert();
		}
		catch (...)
		{
			decomp_flag = true;
			break;
		}
		TMatrixD final_matrix_inverted = final_matrix.Invert();
    for (int l=0; l<9; l++)
    {
      Corrections[l]=0;
      for (int m=0; m<9; m++)
				Corrections[l] += final_matrix_inverted(l,m)*const_matrix[m];
    }

    Combination.pts_with_corrections(Corrections,correctedPt);
    final_chi2 = Combination_unchanged.started_pts(correctedPt) + Combination_unchanged.mass_difference(correctedPt,s_lamda+Corrections[8]);		  
    dif_chi = fabs(final_chi2 - prev_chi);

    float b_left_mass, b_right_mass, left_mass, right_mass;
    Combination.Quartet_final_masses(left_mass, right_mass, Corrections);
    Combination_unchanged.Quartet_final_masses(b_left_mass, b_right_mass, zeroCorrections);
    initialColoronMass[0] = b_left_mass;
    initialColoronMass[1] = b_right_mass;
    correctedColoronMass[0] = left_mass;
    correctedColoronMass[1] = right_mass;


  }while( (iter < nIterations) && (dif_chi > Dif_Chi2) && (final_chi2 >0) && (!TMath::IsNaN(final_chi2))  );
  for(int ipt=0; ipt<8; ipt++)
  {
	  if (correctedPt[ipt] < 0)
		  decomp_flag = true;
  }

  if ( decomp_flag || final_chi2 < 0 || TMath::IsNaN(final_chi2) || iter >= nIterations )
    return -1;
  else
    return 0;
  
}


void VKinematicFit_quartets::setIterations(int iters)
{
  nIterations = iters;
}

void VKinematicFit_quartets::setChi2Dif(double dif)
{
  Dif_Chi2 = dif;
}

void VKinematicFit_quartets::getCorrectedPt(float *new_pt)
{
  for(int i=0; i<8; i++)
    new_pt[i] = correctedPt[i];
}

void VKinematicFit_quartets::getMass_from_Selection(float *masses)
{
  masses[0] = initialColoronMass[0];
  masses[1] = initialColoronMass[1];
}

void VKinematicFit_quartets::getMass_after_Correction(float *masses)
{
  masses[0] = correctedColoronMass[0];
  masses[1] = correctedColoronMass[1];
}

float VKinematicFit_quartets::getChi2Value()
{
  return final_chi2;
}
