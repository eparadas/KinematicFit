#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>



#include "./VKinematicFit_quartets.cpp"
#include "./VKinematicFit_doublets.cpp"


using namespace std;


void v_combs_d(int row, int *hp_combs)
{
	int combs[3][4]={{1,2,3,4},
		{1,3,2,4},
		{1,4,2,3}};
		for(int i=0; i<4; i++)
			hp_combs[i] = combs[row][i]-1;
}




void VKinematicFit_general(int njets, float *o_pt, float *o_eta, float *o_phi, float *o_mass, int *best_combination, float &chi_square_q_, float &chi_square_d_, int *Converged, float *start_pt_, float *intmd_pt_, float *final_pt_, float *f_eta_, float *f_phi_, float *q_ini_mass, float *q_int_mass, float *q_mass, float *d_ini_mass, float *d_int_mass, float *d_mass)
{
	float min_chi_q, min_chi_d;
	float f_f_mass[8], f_mass_[8];
	
	int best_combination_temp[8];
	float f_eta_temp[8], f_phi_temp[8], f_mass_temp[8], start_pt_temp[8], final_pt_temp[8];


	float f_pt[8], f_eta[8], f_phi[8], f_mass[8];
	int temp_hp_comb[8];
  
		
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
		float final_chi;
		int evt_size = 8;
		int fit_q_flag = -1, fit_d_flag = -1;
		for (int ch=0; ch < njets; ch++)
		{
			if ( (o_pt)[ch] < 30 )
			{
				evt_size = ch;
				break;
			}
			else
				evt_size = njets;
		}//for ch

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
	      
// 										cout << "\tProccessing Quartets combination:\t";
										for(int k=0; k<8; k++)
										{
											int k_combi = final_combi[k];
											f_pt[k]=(o_pt)[k_combi];
											start_pts[k]=(o_pt)[k_combi];
											f_eta[k]=(o_eta)[k_combi];
											f_phi[k]=(o_phi)[k_combi];
											f_mass[k]=(o_mass)[k_combi];
// 											cout << k_combi << "\t" ;
										}
// 										cout << "\r" << flush;
	      
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
													f_mass_temp[ji] = f_mass[ji];
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

		if (fit_q_flag == 0)
		{
			min_chi_d = 100000000000000;
	// 		cout << "Performing KinematicFit on doublets...." << endl;
			for(int ihp=0; ihp<3; ihp++)
			{
				v_combs_d(ihp, temp_hp_comb);
				for(int jhp=0; jhp<3; jhp++)
				{
					int temp_jhp_comb[4];
					v_combs_d(jhp, temp_jhp_comb);
					for(int tt=0; tt<4; tt++)
						temp_hp_comb[tt+4] = temp_jhp_comb[tt] + 4;
	
	// 				cout << "Proccessing combination:\t";
					for(int k=0; k<8; k++)
					{
						int k_combi = temp_hp_comb[k];
						f_pt[k]= final_pt_temp[k_combi];
						intmd_pts[k] = final_pt_temp[k_combi];
						start_pts[k] = start_pt_temp[k_combi];
			// 		final_pts[k] = (o_pt)[k_combi];
						f_eta[k]=f_eta_temp[k_combi];
						f_phi[k]=f_phi_temp[k_combi];
						f_f_mass[k]=f_mass_temp[k_combi];
	// 					cout << k_combi << "\t";// << f_pt[k] << "\t" << f_eta[k] << "\t" << f_phi[k] << endl;
					}
	// 				cout << "\r" << flush;
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
								d_ini_mass[cdt] = (temp_jet_ini1 + temp_jet_ini2).M();
								hp_mass_int[cdt] = sqrt( 2 * intmd_pt_[jh] * intmd_pt_[jh+1] * ( cosh(f_eta_[jh] - f_eta_[jh+1]) - cos( f_phi_[jh] - f_phi_[jh+1]) ) );
								cdt++;
							}
	
							for(int jh=0; jh<4; jh++)
								d_int_mass[jh] = hp_mass_int[jh];
			
						}
					}//if fit_status
				}//jhp
			}//ihp
			if ( fit_d_flag != 0 )
			{
	// 			cout << endl << " Because 3C fit didn't return anything, we 're trying M.M.S. " << endl;
				float min_df1 = 100000000000000;
				float min_df2 = 100000000000000;
				for(int ihp=0; ihp<3; ihp++)
				{
					v_combs_d(ihp, temp_hp_comb);
					for(int jhp=0; jhp<3; jhp++)
					{
				
	
						int temp_jhp_comb[4];
						v_combs_d(jhp, temp_jhp_comb);
						for(int tt=0; tt<4; tt++)
							temp_hp_comb[tt+4] = temp_jhp_comb[tt] + 4;
	
	// 					cout << "Proccessing combination:\t";
						for(int k=0; k<8; k++)
						{
							int k_combi = temp_hp_comb[k];
							f_pt[k]= final_pt_temp[k_combi];
							intmd_pts[k] = final_pt_temp[k_combi];
							start_pts[k] = start_pt_temp[k_combi];
				// 		final_pts[k] = (o_pt)[k_combi];
							f_eta_[k]=f_eta_temp[k_combi];
							f_phi_[k]=f_phi_temp[k_combi];
							f_mass_[k]=f_f_mass[k];
	// 						cout << k_combi << "\t" ;
						}
	// 					cout << "\r" << flush;
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
								start_pt_[ji] = start_pts[temp_hp_comb[ji]];
								intmd_pt_[ji] = intmd_pts[temp_hp_comb[ji]];
								final_pt_[ji] = intmd_pt_[temp_hp_comb[ji]];
								f_eta_[ji] = f_eta[temp_hp_comb[ji]];
								f_phi_[ji] = f_phi[temp_hp_comb[ji]];
								f_mass_[ji] = f_f_mass[temp_hp_comb[ji]];
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
		
			chi_square_q_ = min_chi_q;
			chi_square_d_ = min_chi_d;	
		}
		else
			fit_d_flag = fit_q_flag;
}
