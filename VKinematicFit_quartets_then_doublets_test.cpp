#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>


#include "./flat_Class.C"

#include "./VKinematicFit_quartets.cpp"
#include "./VKinematicFit_doublets.cpp"


using namespace std;


void combs_d(int row, int *hp_combs)
{
        int combs[3][4]={{1,2,3,4},
                         {1,3,2,4},
                         {1,4,2,3}};
        for(int i=0; i<4; i++)
                hp_combs[i] = combs[row][i]-1;
}




int main(int argc, char* argv[])
{
    TApplication* rootapp = new TApplication("example",&argc, argv);
    char outputroot[64], outputtxt[64];
    int num_file = atoi(argv[3]);
    sprintf(outputroot, "vkinematicfit_%d.root", num_file);
    sprintf(outputtxt, "vkinematicfit_%d.txt", num_file);
    TFile *outf = new TFile(outputroot,"recreate");
    
    TTree *outTree = new TTree("outTree","events");
    
  char tm;
//      ofstream out("chi_values.txt");
  ofstream out_res(outputtxt);
  flat_Class t1;
  float left_mass, right_mass, min_chi_q, min_chi_d, best_min_chi ;
  float m4jave, m2jave, f_f_mass[8], f_mass_[8], q_ini_mass[2], q_int_mass[2], q_mass[2], chi_square_q_, chi_square_d_, q_average_mass, d_ini_mass[4], d_int_mass[4], d_mass[4], d1_average_mass, d2_average_mass;

  float Corrections[9];
  int best_combination[8], best_combination_temp[8], i, Converged[2];
  float pt_fd[8], eta_fd[8], phi_fd[8],start_pt_[8], intmd_pt_[8], final_pt_[8], f_eta_[8], f_phi_[8], f_eta_temp[8], f_phi_temp[8], start_pt_temp[8], final_pt_temp[8];
  long int count_combs = 0, count_ninv = 0; 
  outTree->Branch("event", &i, "i/I");
  outTree->Branch("best_comb", &best_combination, "i[8]/I");
  outTree->Branch("chi_square_quartets", &chi_square_q_, "chi_square_q_/F");
  outTree->Branch("chi_square_doublets", &chi_square_d_, "chi_square_d_/F");
  outTree->Branch("start_pt", &start_pt_, "start_pt_[8]/F");
  outTree->Branch("intmd_pt", &intmd_pt_, "intmd_pt_[8]/F");
  outTree->Branch("final_pt", &final_pt_, "final_pt_[8]/F");
  outTree->Branch("f_eta", &f_eta_, "f_eta_[8]/F");
  outTree->Branch("f_phi", &f_phi_, "f_phi_[8]/F");
  outTree->Branch("quartet_ini_mass", &q_ini_mass, "q_ini_mass[2]/F");
  outTree->Branch("quartet_int_mass", &q_int_mass, "q_int_mass[2]/F");
  outTree->Branch("quartet_mass", &q_mass, "q_mass[2]/F");
  outTree->Branch("doublet_ini_mass", &d_ini_mass, "d_ini_mass[4]/F");
  outTree->Branch("doublet_int_mass", &d_int_mass, "d_int_mass[4]/F");
  outTree->Branch("doublet_mass", &d_mass, "d_mass[4]/F");
  


  int start_for = atoi(argv[1]);
  int end_for = atoi(argv[2]);

  float f_pt[8], f_eta[8], f_phi[8], sigma[8], f_mass[8];
  int temp_hp_comb[8];
  
  for (i=start_for; i<end_for; i++)
  {
  // 	cout <<"Event " << i << endl;
  //	out<<endl;
  t1.GetEntry(i);
  if ( t1.total_ht < 1000 || (*t1.pt)[7] < 30 )
  continue;

  cout << "Event: " << i << "\nPerforming KinematicFit on quartets...." << endl;
//----initialize output variables--------------------  
  for(int jk=0; jk<8; jk++)
  {
    start_pt_[jk]	= -999;
    intmd_pt_[jk] 	= -999;
    final_pt_[jk] 	= -999;
    f_eta_[jk]    	= -999;
    f_phi_[jk]	 	= -999;
    if ( jk < 2 )
    {
      q_ini_mass[jk]	= -999;
      q_int_mass[jk]	= -999;
      Converged[jk] 	= -999;
      q_mass[jk]    	= -999;
    }
    if (jk < 4)
    {
      d_ini_mass[jk] 	= -999;
      d_int_mass[jk]	= -999;
      d_mass[jk]	= -999;
    }
    
    chi_square_d_	= -999;
    chi_square_q_ 	= -999;
  }
//--------------------------------------------------  
  min_chi_q = 100000000000000;
  best_min_chi = 10000;
  float final_chi;
  int evt_size;
  int fit_q_flag = -1, fit_d_flag = -1;
  for (int ch=0; ch < (*t1.pt).size(); ch++)
  {
  if ( (*t1.pt)[ch] < 30 )
  {
    evt_size = ch;
    break;
  }
  else
    evt_size = (*t1.pt).size();
  }//for ch

  if ( evt_size < 8 )
	continue;
  float start_pts[8], final_pts[8], intmd_pts[8];
  int final_combi[8];
  bool flag;
  vector<int> nums, temp_num;

  nums.clear();
  for(int lii=0; lii<evt_size; lii++)
  nums.push_back(lii);

  for (int j1=0; j1<evt_size-7; j1++)
  {
  for (int j2=j1+1; j2<evt_size-2; j2++)
  {
  if ( j2 == j1 )
    continue;

  for (int j3=j2+1; j3<evt_size-1; j3++)
  {
    if ( j3 == j2 || j3 == j1 )
      continue;
    for (int j4=j3+1; j4<evt_size; j4++)
    {
      if ( j4 == j3 || j4 == j2 || j4 == j1 )
	continue;
      
      nums.clear();
      for(int lii=j1; lii<evt_size; lii++)
	nums.push_back(lii);
      
  // 		cout << j1 << "\t" << j2 << "\t" << j3 << "\t" << j4 << "\t";
  // 	      out << j1 << "\t" << j2 << "\t" << j3 << "\t" << j4 << "\t";
      final_combi[0] = j1;
      final_combi[1] = j2;
      final_combi[2] = j3;
      final_combi[3] = j4;
      
      int count2=0;
      for (int k=0; k<evt_size; k++)
      {
	flag = false;
	for (int j=0; j<4; j++)
	{
	  if ( final_combi[j] == k || count2 == 4 )
	    flag = true;
	}
	
	if ( flag )
	  continue;
  // 		  cout << k << "\t";
  // 		out << k << "\t";
  // 		if ( k < final_combi[count2] )
  // 		  continue;
      }
  //--------erase first selected numbers-------------
      for(int jk=0; jk<4; jk++)
      {
	temp_num.clear();
	temp_num.push_back(final_combi[jk]);
	
	vector<int>::iterator it;
	it = search(nums.begin(), nums.end(), temp_num.begin(), temp_num.end());
	if (it!=nums.end())
	{
	  int pos(it-nums.begin());
	  nums.erase(nums.begin()+pos);
	}
      }
  //---------------------------------------------------	      
      for (int j5=0; j5<evt_size-4 -j1; j5++)
      {
	for (int j6 = j5+1; j6<evt_size-4 -j1; j6++)
	{
	  if ( j6 == j5 )
	    continue;
	  
	  for (int j7 = j6+1; j7<evt_size-4 -j1; j7++)
	  {
	    if ( j7==j6 || j7 == j5 )
	      continue;
	    
	    for (int j8 = j7+1; j8<evt_size-4 -j1; j8++)
	    {
	      if ( j8 == j7 || j8 == j6 || j8 == j5 )
		continue;
	      
	      final_combi[4] = nums[j5];
	      final_combi[5] = nums[j6];
	      final_combi[6] = nums[j7];
	      final_combi[7] = nums[j8];
	      
	      int lead_count(0);
	      for(int jk=0; jk<8; jk++)
	      {
		if ( final_combi[jk] < 5 )
		  lead_count++;
	      }
	      if ( lead_count < 5 )
		continue;
	      
	      cout << "\tProccessing Quartets combination:\t";
	      for(int k=0; k<8; k++)
	      {
		    int k_combi = final_combi[k];
		    f_pt[k]=(*t1.pt)[k_combi];
		    start_pts[k]=(*t1.pt)[k_combi];
		    f_eta[k]=(*t1.eta)[k_combi];
		    f_phi[k]=(*t1.phi)[k_combi];
		    f_mass[k]=(*t1.mass)[k_combi];
		    cout << k_combi << "\t" ;
	      }
	      cout << "\r" << flush;
	      
	      VKinematicFit_quartets Combination(f_pt, f_eta, f_phi);
	      int fit_status_q = Combination.ApplyFit();
	      
	      if ( fit_status_q == 0 )
	      {
		fit_q_flag = fit_status_q;
		Converged[0] = 0;
		final_chi = Combination.getChi2Value();
		if ( final_chi < min_chi_q )
		{
		  min_chi_q = final_chi;
		
		  Combination.getCorrectedPt(final_pts);
		  float in_masses[2], fi_masses[2];
		  Combination.getMass_from_Selection(in_masses);
		  Combination.getMass_after_Correction(fi_masses);
		  q_int_mass[0] = fi_masses[0];
		  q_int_mass[1] = fi_masses[1];
		  
		  for(int ji = 0; ji < 8; ji++)
		  {
		    best_combination_temp[ji] = final_combi[ji];
		    start_pt_temp[ji] = start_pts[ji];
		    final_pt_temp[ji] = final_pts[ji];
		    f_eta_temp[ji] = f_eta[ji];
		    f_phi_temp[ji] = f_phi[ji];
		  }
		  
		  TLorentzVector col_1(0,0,0,0);
		  TLorentzVector col_2(0,0,0,0);
		  for(int jc=0; jc<4; jc++)
		  {
		    TLorentzVector temp_jet(0,0,0,0);
		    temp_jet.SetPtEtaPhiM(start_pts[jc], f_eta[jc], f_phi[jc], f_mass[jc]);
		    col_1 += temp_jet;
		  }
		  for(int jc=4; jc<8; jc++)
		  {
		    TLorentzVector temp_jet(0,0,0,0);
		    temp_jet.SetPtEtaPhiM(start_pts[jc], f_eta[jc], f_phi[jc], f_mass[jc]);
		    col_2 += temp_jet;
		  }
		  
		  q_ini_mass[0] = col_1.M();
		  q_ini_mass[1] = col_2.M();
		}
	      }//if fit_status
	      
	      }//j8
	  }//j7
	}//j6
      }//j5
      }//j4
  // 	cout << endl;
    }//j3
  }//j2
  }//j1

  min_chi_d = 100000000000000;
  cout << "Performing KinematicFit on doublets...." << endl;
  for(int ihp=0; ihp<3; ihp++)
  {
    combs_d(ihp, temp_hp_comb);
    for(int jhp=0; jhp<3; jhp++)
    {
		
    // 			float f_pt[8], f_eta[8], f_phi[8], sigma[8];
      int temp_jhp_comb[4];
      combs_d(jhp, temp_jhp_comb);
      for(int tt=0; tt<4; tt++)
	temp_hp_comb[tt+4] = temp_jhp_comb[tt] + 4;

      cout << "Proccessing combination:\t";
      for(int k=0; k<8; k++)
      {
	  int k_combi = temp_hp_comb[k];
	  f_pt[k]= final_pt_temp[k_combi];
	  intmd_pts[k] = final_pt_temp[k_combi];
	  start_pts[k] = start_pt_temp[k_combi];
    // 		final_pts[k] = (*t1.pt)[k_combi];
	  f_eta[k]=f_eta_temp[k_combi];
	  f_phi[k]=f_phi_temp[k_combi];
	  cout << k_combi << "\t";// << f_pt[k] << "\t" << f_eta[k] << "\t" << f_phi[k] << endl;
      }
      cout << "\r" << flush;
//       cout << endl;
      VKinematicFit_doublets Combination(f_pt, f_eta, f_phi);
      int fit_status_d = Combination.ApplyFit();
      
      if ( fit_status_d == 0 )
      {
	fit_d_flag = fit_status_d;
	Converged[1] = 0;
	final_chi = Combination.getChi2Value();
	if ( final_chi < min_chi_d )
	{
	  min_chi_d = final_chi;
	
	  Combination.getCorrectedPt(final_pts);
	  float in_masses_1[2], in_masses_2[2], fi_masses_1[2], fi_masses_2[2], q_masses[2];
	  Combination.getMass_from_Selection(in_masses_1, in_masses_2);
	  Combination.getMass_after_Correction(fi_masses_1,fi_masses_2);
	  Combination.getMass_for_Quartets(q_masses);
	  
	  d_mass[0] = fi_masses_1[0];
	  d_mass[1] = fi_masses_1[1];
	  d_mass[2] = fi_masses_1[0];
	  d_mass[3] = fi_masses_1[1];
	  q_mass[0] = q_masses[0];
	  q_mass[1] = q_masses[1];
	  
	  
	  for(int ji = 0; ji < 8; ji++)
	  {
	    best_combination[ji] = best_combination_temp[temp_hp_comb[ji]];
	    start_pt_[ji] = start_pts[ji];
	    intmd_pt_[ji] = intmd_pts[ji];
	    final_pt_[ji] = final_pts[ji];
	    f_eta_[ji] = f_eta[ji];
	    f_phi_[ji] = f_phi[ji];
	    f_mass_[ji] = f_f_mass[ji];
	  }
		
		int cdt = 0;
		float hp_mass_int[4];
		TLorentzVector hp_ini[4];
		for(int jh=0; jh<8; jh=jh+2)
		{
			TLorentzVector temp_jet_ini1(0,0,0,0);
			TLorentzVector temp_jet_ini2(0,0,0,0);

			temp_jet_ini1.SetPtEtaPhiM(start_pt_[jh+1], f_eta_[jh+1], f_phi_[jh+1], f_mass_[jh+1]);
			temp_jet_ini2.SetPtEtaPhiM(start_pt_[jh], f_eta_[jh], f_phi_[jh], f_mass_[jh]);
			hp_ini[cdt] = (temp_jet_ini1 + temp_jet_ini2);
			
			hp_mass_int[cdt] = sqrt( 2 * intmd_pt_[jh] * intmd_pt_[jh+1] * ( cosh(f_eta_[jh] - f_eta_[jh+1]) - cos( f_phi_[jh] - f_phi_[jh+1]) ) );
			cdt++;
		}

		for(int jh=0; jh<4; jh++)
		{
			d_ini_mass[jh] = hp_ini[jh].M();
			d_int_mass[jh] = hp_mass_int[jh];
		}
		
	}
      }//if fit_status
    }//jhp
  }//ihp
  if ( fit_d_flag != 0 )
  {
    cout << endl << " Because 3C fit didn't return anything, we 're trying M.M.S. " << endl;
    float min_df1 = 100000000000000;
    float min_df2 = 100000000000000;
    float min_df3 = 100000000000000;
    for(int ihp=0; ihp<3; ihp++)
    {
      combs_d(ihp, temp_hp_comb);
      for(int jhp=0; jhp<3; jhp++)
      {
		  
      // 			float f_pt[8], f_eta[8], f_phi[8], sigma[8];
	int temp_jhp_comb[4];
	combs_d(jhp, temp_jhp_comb);
	for(int tt=0; tt<4; tt++)
	  temp_hp_comb[tt+4] = temp_jhp_comb[tt] + 4;

	cout << "Proccessing combination:\t";
	for(int k=0; k<8; k++)
	{
	    int k_combi = temp_hp_comb[k];
	    f_pt[k]= final_pt_temp[k_combi];
	    intmd_pts[k] = final_pt_temp[k_combi];
	    start_pts[k] = start_pt_temp[k_combi];
      // 		final_pts[k] = (*t1.pt)[k_combi];
	    f_eta_[k]=f_eta_temp[k_combi];
	    f_phi_[k]=f_phi_temp[k_combi];
	    f_mass_[k]=f_f_mass[k];
	    cout << k_combi << "\t" ;
	}
	cout << "\r" << flush;
	int cdt = 0;
	TLorentzVector hp[4];
	for(int jh=0; jh<8; jh=jh+2)
	{
	  TLorentzVector temp_jet(0,0,0,0);
	  temp_jet.SetPtEtaPhiM(f_pt[jh+1], f_eta_[jh+1], f_phi_[jh+1], f_mass_[jh+1]);
	  hp[cdt].SetPtEtaPhiM(f_pt[jh], f_eta_[jh], f_phi_[jh], f_mass_[jh]);
	  hp[cdt] += temp_jet; 
	}
	
	float temp_min_dif1 = fabs(hp[0].M() - hp[1].M());
	float temp_min_dif2 = fabs(hp[2].M() - hp[3].M());
	
	if ( temp_min_dif1 < min_df1 )
	{
	  min_df1 = temp_min_dif1;
	  d_mass[0] = hp[0].M();
	  d_mass[1] = hp[1].M();
	  
	  for(int ji = 0; ji < 4; ji++)
	  {
	    best_combination[ji] = best_combination_temp[temp_hp_comb[ji]];
	    start_pt_[ji] = start_pts[ji];
	    intmd_pt_[ji] = intmd_pts[ji];
	    final_pt_[ji] = intmd_pt_[ji];
	    f_eta_[ji] = f_eta[ji];
	    f_phi_[ji] = f_phi[ji];
	    f_mass_[ji] = f_f_mass[ji];
	  }
	}
	
	if ( temp_min_dif2 < min_df2 )
	{
	  min_df2 = temp_min_dif2;
	  d_mass[2] = hp[2].M();
	  d_mass[3] = hp[3].M();
	  
	  for(int ji = 4; ji < 8; ji++)
	  {
	    best_combination[ji] = best_combination_temp[temp_hp_comb[ji]];
	    start_pt_[ji] = start_pts[ji];
	    intmd_pt_[ji] = intmd_pts[ji];
	    final_pt_[ji] = intmd_pt_[ji];
	    f_eta_[ji] = f_eta[ji];
	    f_phi_[ji] = f_phi[ji];
	    f_mass_[ji] = f_f_mass[ji];
	  }
	}
      }//jhp
    }//ihp
    q_mass[0] = q_ini_mass[0];
    q_mass[1] = q_ini_mass[1];
  }// if !fit_d_flag
	
	
    for(int tt=0; tt<8; tt++)
    {
      cout << endl << setprecision(10) << endl << best_combination[tt] << "\t" << start_pt_[tt] << "\t" <<intmd_pt_[tt] << "\t" << final_pt_[tt] << "\t" << f_eta_[tt] << "\t" << f_phi_[tt] << endl;
    }
    cout << endl<<endl;
    for(int ti=0; ti<(*t1.pt).size(); ti++)
      cout << ti << "\t"<< (*t1.pt)[ti] << "\t" << (*t1.eta)[ti] << "\t" << (*t1.phi)[ti] << endl;
    cout << "Masses:\t" << q_mass[0] << "\t" << q_mass[0] << "\t" << d_mass[0] << "\t" << d_mass[1] << "\t" << d_mass[2] << "\t" << d_mass[3] << endl;
    cout << "Fit status:\t" << fit_q_flag << "\t" << fit_d_flag << endl;
    
    out_res << "Event: " << i << endl << "\tCombination is: ";
    for (int cj=0; cj<8; cj++)
      out_res << best_combination[cj] << "\t";
    chi_square_q_ = min_chi_q;
    chi_square_d_ = min_chi_d;
//     q_average_mass = (q_left_mass + q_right_mass)/2.;
    out_res << endl << "with min_chi = " << min_chi_q << "\t" << min_chi_d << endl << " m4j_left = " << q_mass[0] << " m4j_right = " << q_mass[1] << endl << "m2jave " << d_mass[0] << "\t" << d_mass[1] << "\t" << d_mass[2] << "\t" << d_mass[3] << endl;
    out_res << endl;
    cout << "chi_values:\t" << chi_square_q_ << "\t" << chi_square_d_ << endl;
  
    outTree->Fill();

//   cin>>tm;

  }// for i

  outf -> cd();
  outTree -> Write();
  outf -> Close(); 
  cout << "Finally non inverted matrices are: " << count_ninv << " and total number of combinations is: " << count_combs << endl;
  rootapp->Run();
      
}