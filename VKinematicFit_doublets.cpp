#include "./VKinematicFit_doublets.h"
// #include "./nVDecompLU.cpp"
#include "./KinFit_class_with_constraints_quadrupal_only_doublets.cpp"

VKinematicFit_doublets::VKinematicFit_doublets(float *f_pt, float *f_eta, float *f_phi)
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

void VKinematicFit_doublets::initialize()
{
  for(int k=0; k<3; k++)
    s_lamda[k] = 1.0;
  
  nIterations = 200;
  countedIterations = 0;
  Dif_Chi2 = 1e-12;
  final_chi2 = 100000000;
}

int VKinematicFit_doublets::ApplyFit()
{
  initialize();
  bool decomp_flag = false;

  Pts_dd Combination_unchanged( initialPt, initialEta, initialPhi, s_lamda );

    
  int iter(0);
  double prev_chi, dif_chi;
  do{
    iter++;
    prev_chi = final_chi2;
    for (int k=0; k<8; k++)	
      initialPt[k] += Corrections[k];
      
    for(int k=0; k<3; k++)
      s_lamda[k] += Corrections[k+8];

    
    Pts_dd Combination( initialPt, initialEta, initialPhi, s_lamda );
    
    double const_matrix[11];
    for (int l=0; l<11; l++)
      const_matrix[l] = (Combination.Derived_pt(l,11) + Combination.Derived_1st_doublet_term(l,11) + Combination.Derived_3rd_doublet_term(l,11) + Combination.Derived_d1_lamda_term(l,11) + Combination.Derived_d2_lamda_term(l,11) + Combination.Derived_d3_lamda_term(l,11));


    TMatrixD final_matrix(11,11);
    for (int l=0; l<11; l++)
    {
      for (int m=0; m<11; m++)
	final_matrix(l,m) = (Combination.Derived_pt(l,m) + Combination.Derived_1st_doublet_term(l,m) + Combination.Derived_3rd_doublet_term(l,m) + Combination.Derived_d1_lamda_term(l,m) + Combination.Derived_d2_lamda_term(l,m) + Combination.Derived_d3_lamda_term(l,m));
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
		
    for (int l=0; l<11; l++)
    {
      Corrections[l]=0;
      for (int m=0; m<11; m++)
				Corrections[l] += final_matrix_inverted(l,m)*const_matrix[m];
    }

    Combination.pts_with_corrections(Corrections, correctedPt);
    final_chi2 = Combination_unchanged.started_pts(correctedPt) + Combination_unchanged.mass_difference1(correctedPt,s_lamda[0]+Corrections[8]) + Combination_unchanged.mass_difference2(correctedPt,s_lamda[1]+Corrections[9]) + Combination_unchanged.mass_difference3(correctedPt,s_lamda[2]+Corrections[10]);
    dif_chi = fabs(final_chi2 - prev_chi);

    float b_left_mass_1, b_left_mass_2, b_right_mass_1, b_right_mass_2, left_mass_1, left_mass_2, right_mass_1, right_mass_2, q_left_mass, q_right_mass;
    Combination.Doublet1_final_masses(left_mass_1, right_mass_1, Corrections);
    Combination.Doublet2_final_masses(left_mass_2, right_mass_2, Corrections);
    Combination_unchanged.Doublet1_final_masses(b_left_mass_1, b_right_mass_1, zeroCorrections);
    Combination_unchanged.Doublet2_final_masses(b_left_mass_2, b_right_mass_2, zeroCorrections);
    Combination.Quartet_final_masses(q_left_mass, q_right_mass, Corrections);
    
    initialDoublet_1Mass[0] = b_left_mass_1;
    initialDoublet_1Mass[1] = b_right_mass_1;
    initialDoublet_2Mass[0] = b_left_mass_2;
    initialDoublet_2Mass[1] = b_right_mass_2;
    correctedDoublet_1Mass[0] = left_mass_1;
    correctedDoublet_1Mass[1] = right_mass_1;
    correctedDoublet_2Mass[0] = left_mass_2;
    correctedDoublet_2Mass[1] = right_mass_2;
    correctedQuartetmass[0] = q_left_mass;
    correctedQuartetmass[1] = q_right_mass;
    
  }while( (iter < nIterations) && (dif_chi > Dif_Chi2) && (final_chi2 >0) && (!TMath::IsNaN(final_chi2))  );
  for(int ipt=0; ipt<8; ipt++)
  {
	  if (correctedPt[ipt] < 0)
		  decomp_flag = true;
  }
  

  if ( decomp_flag || final_chi2 < 0 || TMath::IsNaN(final_chi2) || iter >= nIterations)
    return -1;
  else
    return 0;

}


void VKinematicFit_doublets::setIterations(int iters)
{
  nIterations = iters;
}

void VKinematicFit_doublets::setChi2Dif(double dif)
{
  Dif_Chi2 = dif;
}

void VKinematicFit_doublets::getCorrectedPt(float *new_pt)
{
  for(int i=0; i<8; i++)
    new_pt[i] = correctedPt[i];
}

void VKinematicFit_doublets::getMass_from_Selection(float *masses_1, float *masses_2)
{
  masses_1[0] = initialDoublet_1Mass[0];
  masses_1[1] = initialDoublet_1Mass[1];
  masses_2[0] = initialDoublet_2Mass[0];
  masses_2[1] = initialDoublet_2Mass[1];
  
}

void VKinematicFit_doublets::getMass_after_Correction(float *masses_1, float *masses_2)
{
  masses_1[0] = correctedDoublet_1Mass[0];
  masses_1[1] = correctedDoublet_1Mass[1];
  masses_2[0] = correctedDoublet_2Mass[0];
  masses_2[1] = correctedDoublet_2Mass[1];
  
}

void VKinematicFit_doublets::getMass_for_Quartets(float *masses)
{
  masses[0] = correctedQuartetmass[0];
  masses[1] = correctedQuartetmass[1];
}

float VKinematicFit_doublets::getChi2Value()
{
  return final_chi2;
}
