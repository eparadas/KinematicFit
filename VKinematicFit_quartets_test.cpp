#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "./flat_Class.C"

#include "./VKinematicFit_quartets.cpp"


using namespace std;

int main(int argc, char* argv[])
{
    TApplication* rootapp = new TApplication("example",&argc, argv);

    TFile *outf = new TFile("vkinematicfit_test.root","recreate");
    
    TTree *outTree = new TTree("outTree","events");
    
  char tm;
//      ofstream out("chi_values.txt");
      ofstream out_res("vkinematicfit_test.txt");
      flat_Class t1;
      float left_mass, right_mass, min_chi, best_min_chi ;
      float m4jave, m2jave, q_left_mass, q_right_mass, chi_square_, q_average_mass, d1_left_mass, d1_right_mass, d2_left_mass, d2_right_mass, d1_average_mass, d2_average_mass;

       float Corrections[9];
       int best_combination[8],i;
       float pt_fd[8], eta_fd[8], phi_fd[8],start_pt_[8], final_pt_[8], f_eta_[8], f_phi_[8];
       long int count_combs = 0, count_ninv = 0; 
       outTree->Branch("event", &i, "i/I");
       outTree->Branch("best_comb", &best_combination, "i[8]/I");
       outTree->Branch("chi_square", &chi_square_, "chi_square_/F");
       outTree->Branch("start_pt", &start_pt_, "start_pt[8]_/F");
       outTree->Branch("final_pt", &final_pt_, "final_pt_[8]/F");
       outTree->Branch("f_eta", &f_eta_, "f_eta_[8]/F");
       outTree->Branch("f_phi", &f_phi_, "f_phi_[8]/F");
       outTree->Branch("quartet_lmass", &q_left_mass, "q_left_mass/F");
       outTree->Branch("quartet_rmass", &q_right_mass, "q_right_mass/F");
       outTree->Branch("quartet_mass_ave", &q_average_mass, "q_average_mass/F");
       outTree->Branch("doublet1_lmass", &d1_left_mass, "d1_left_mass/F");
       outTree->Branch("doublet1_rmass", &d1_right_mass, "d1_right_mass/F");
       outTree->Branch("doublet1_mass_ave", &d1_average_mass, "d1_average_mass/F");
       outTree->Branch("doublet2_lmass", &d2_left_mass, "d2_left_mass/F");
       outTree->Branch("doublet2_rmass", &d2_right_mass, "d2_right_mass/F");
       outTree->Branch("doublet2_mass_ave", &d2_average_mass, "d2_average_mass/F");
       
       int start_for = atoi(argv[1]);
       int end_for = atoi(argv[2]);

       float f_pt[8], f_eta[8], f_phi[8], sigma[8];
      
      for (i=start_for; i<end_for; i++)
      {
// 	cout <<"Event " << i << endl;
//	out<<endl;
	t1.GetEntry(i);
 	if ( t1.total_ht < 1000 || (*t1.pt)[7] < 30 )
 	  continue;
	
	min_chi = 100000000000000;
	best_min_chi = 10000;
	float final_chi = 0;
	int evt_size;
	for (int ch=0; ch < (*t1.pt).size(); ch++)
	{
	  if ( (*t1.pt)[ch] < 30 )
	  {
	    evt_size = ch+1;
	    break;
	  }
	  else
	    evt_size = (*t1.pt).size();
	}//for ch
	
	if ( evt_size < 8 )
		continue;
	float start_pts[8], final_pts[8];
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
		      
		      cout << "Event: " << i << "\tProccessing Quartets combination:\t";
		      for(int k=0; k<8; k++)
		      {
			    int k_combi = final_combi[k];
			    f_pt[k]=(*t1.pt)[k_combi];
			    start_pts[k]=(*t1.pt)[k_combi];
			    f_eta[k]=(*t1.eta)[k_combi];
			    f_phi[k]=(*t1.phi)[k_combi];
			    cout << k_combi << "\t" ;
		      }
		      cout << "\r" << flush;
		      
		      VKinematicFit_quartets Combination(f_pt, f_eta, f_phi);
		      int fit_status = Combination.ApplyFit();
		      
		      if ( fit_status == 0 )
		      {
			final_chi = Combination.getChi2Value();
			if ( final_chi < min_chi )
			{
			  min_chi = final_chi;
			
			  Combination.getCorrectedPt(final_pts);
			  float in_masses[2], fi_masses[2];
			  Combination.getMass_from_Selection(in_masses);
			  Combination.getMass_after_Correction(fi_masses);
			  q_left_mass = fi_masses[0];
			  q_right_mass= fi_masses[1];
			  
			  for(int ji = 0; ji < 8; ji++)
			  {
			    best_combination[ji] = final_combi[ji];
			    start_pt_[ji] = start_pts[ji];
			    final_pt_[ji] = final_pts[ji];
			    f_eta_[ji] = f_eta[ji];
			    f_phi_[ji] = f_phi[ji];
			  }
			}
		      }
		      
		      }//j8
		  }//j7
		}//j6
	      }//j5
	      }//j4
     // 	cout << endl;
	    }//j3
	  }//j2
	}//j1
	
	for(int tt=0; tt<8; tt++)
	{
	  cout << endl << setprecision(10) << endl << best_combination[tt] << "\t" << f_pt[tt] << "\t" << start_pt_[tt] << "\t" << final_pt_[tt] << "\t" << f_eta_[tt] << "\t" << f_phi_[tt] << endl;
	}
	cout << endl<<endl;
	for(int ti=0; ti<(*t1.pt).size(); ti++)
	  cout << ti << "\t"<< (*t1.pt)[ti] << "\t" << (*t1.eta)[ti] << "\t" << (*t1.phi)[ti] << endl;
	cout << "Masses:\t" << q_left_mass << "\t" << q_right_mass << endl;
	out_res << "Event: " << i << endl << "\tCombination is: ";
	for (int cj=0; cj<8; cj++)
	  out_res << best_combination[cj] << "\t";
	chi_square_ = min_chi;
	q_average_mass = (q_left_mass + q_right_mass)/2.;
	out_res << endl << "with Best min_chi = " << min_chi << endl << " m4j_left = " << q_left_mass << " m4j_right = " << q_right_mass << endl << " m4jave = " << q_average_mass << " and m2jave " << d1_average_mass << "\t" << d2_average_mass << endl;
	out_res << endl;
	outTree->Fill();

//	cin>>tm;

      }// for i
      
      outf -> cd();
      outTree -> Write();
      outf -> Close(); 
	cout << "Finally non inverted matrices are: " << count_ninv << " and total number of combinations is: " << count_combs << endl;
       rootapp->Run();
      
}