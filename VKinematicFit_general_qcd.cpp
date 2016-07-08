#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <TLorentzVector.h>
#include <TApplication.h>


#include "./flat_qcd.C"
#include "./VKinematicFit_general.cpp"

int main(int argc, char* argv[])
{
	TApplication *rootapp = new TApplication("kinematicfit", &argc, argv);
	
	flat_qcd t1;
	
	char outputroot[64], outputtxt[64];
	int num_file = atoi(argv[3]);
	sprintf(outputroot, "vkinematicfit_qcd_%d.root", num_file);
	sprintf(outputtxt, "vkinematicfit_qcd_%d.txt", num_file);
	TFile *outf = new TFile(outputroot,"recreate");
    
	TTree *outTree = new TTree("outTree","events");
	
	char tm;
//      ofstream out("chi_values.txt");
	ofstream out_res(outputtxt);
	
	float left_mass, right_mass, min_chi_q, min_chi_d, best_min_chi ;
	float m4jave, m2jave, f_f_mass[8], f_mass_[8], q_ini_mass_mkf_[2], q_int_mass_mkf_[2], q_mass_mkf_[2], chi_square_q_mkf_, chi_square_d_mkf_, q_average_mass, d_ini_mass_mkf_[4], d_int_mass_mkf_[4], d_mass_mkf_[4], d1_average_mass, d2_average_mass;

	float Corrections[9];
	int i, best_combination_mkf[8], best_combination_mkf_temp[8], evt, Converged_mkf_[2];
	float pt_fd[8], eta_fd[8], phi_fd[8],start_pt_mkf_[8], intmd_pt_mkf_[8], final_pt_mkf_[8], f_eta_mkf_[8], f_phi_mkf_[8], f_eta_mkf_temp[8], f_phi_mkf_temp[8], start_pt_mkf_temp[8], final_pt_mkf_temp[8];
	long int count_combs = 0, count_ninv = 0; 
	
	outTree->Branch("event_mkf", &i, "i/I");
	outTree->Branch("best_comb_mkf", &best_combination_mkf, "&best_combination_mkf[8]/I");
	outTree->Branch("chi_square_quartets_mkf", &chi_square_q_mkf_, "chi_square_q_mkf_/F");
	outTree->Branch("chi_square_doublets_mkf", &chi_square_d_mkf_, "chi_square_d_mkf_/F");
	outTree->Branch("converge_mkf", &Converged_mkf_, "Converged_mkf_[2]/I");
	outTree->Branch("start_pt_mkf", &start_pt_mkf_, "start_pt_mkf_[8]/F");
	outTree->Branch("intmd_pt_mkf", &intmd_pt_mkf_, "intmd_pt_mkf_[8]/F");
	outTree->Branch("final_pt_mkf", &final_pt_mkf_, "final_pt_mkf_[8]/F");
	outTree->Branch("f_eta_mkf", &f_eta_mkf_, "f_eta_mkf_[8]/F");
	outTree->Branch("f_phi_mkf", &f_phi_mkf_, "f_phi_mkf_[8]/F");
	outTree->Branch("quartet_ini_mass_mkf", &q_ini_mass_mkf_, "q_ini_mass_mkf_[2]/F");
	outTree->Branch("quartet_int_mass_mkf", &q_int_mass_mkf_, "q_int_mass_mkf_[2]/F");
	outTree->Branch("quartet_mass_mkf", &q_mass_mkf_, "q_mass_mkf_[2]/F");
	outTree->Branch("doublet_ini_mass_mkf", &d_ini_mass_mkf_, "d_ini_mass_mkf_[4]/F");
	outTree->Branch("doublet_int_mass_mkf", &d_int_mass_mkf_, "d_int_mass_mkf_[4]/F");
	outTree->Branch("doublet_mass_mkf", &d_mass_mkf_, "d_mass_mkf_[4]/F");
  


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

		evt = t1.evtNo;
		cout << "i = " << i << "\tEvent: " << evt << "\nPerforming KinematicFit on quartets...." << endl;
		
		const int njets = (*t1.pt).size();
		float ev_pt[njets], ev_eta[njets], ev_phi[njets], ev_mass[njets];
		
		for (int j=0; j<njets; j++)
		{
			ev_pt[j] = (*t1.pt)[j];
			ev_eta[j] = (*t1.eta)[j];
			ev_phi[j] = (*t1.phi)[j];
			ev_mass[j] = (*t1.mass)[j];
		}
		
		VKinematicFit_general(njets, ev_pt,  ev_eta,  ev_phi,  ev_mass, best_combination_mkf,  chi_square_q_mkf_,  chi_square_d_mkf_,  Converged_mkf_,  start_pt_mkf_,  intmd_pt_mkf_,  final_pt_mkf_,  f_eta_mkf_,  f_phi_mkf_,  q_ini_mass_mkf_,  q_int_mass_mkf_,  q_mass_mkf_,  d_ini_mass_mkf_,  d_int_mass_mkf_,  d_mass_mkf_);
		
		outTree->Fill();
		
//		for(int tt=0; tt<8; tt++)
//		{
//			cout << endl << setprecision(10) << endl << best_combination_mkf[tt] << "\t" << start_pt_mkf_[tt] << "\t" <<intmd_pt_mkf_[tt] << "\t" << final_pt_mkf_[tt] << "\t" << f_eta_mkf_[tt] << "\t" << f_phi_mkf_[tt] << endl;
//		}
//		cout << endl<<endl;
//		for(int ti=0; ti<(*t1.pt).size(); ti++)
//			cout << ti << "\t"<< (*t1.pt)[ti] << "\t" << (*t1.eta)[ti] << "\t" << (*t1.phi)[ti] << endl;
//		cout << "Masses:\t" << q_mass_mkf_[0] << "\t" << q_mass_mkf_[0] << "\t" << d_mass_mkf_[0] << "\t" << d_mass_mkf_[1] << "\t" << d_mass_mkf_[2] << "\t" << d_mass_mkf_[3] << endl;
//		cout << "Fit status:\t" << Converged_mkf_[0] << "\t" << Converged_mkf_[1] << endl;
    
		out_res << "Event: " << i << endl << "\tCombination is: ";
		for (int cj=0; cj<8; cj++)
			out_res << best_combination_mkf[cj] << "\t";
// 		chi_square_q_mkf = min_chi_q;
// 		chi_square_d_mkf_ = min_chi_d;
//     q_average_mass = (q_left_mass + q_right_mass)/2.;
		out_res << endl << "with min_chi = " << chi_square_q_mkf_ << "\t" << chi_square_d_mkf_ << endl << " m4j_left = " << q_mass_mkf_[0] << " m4j_right = " << q_mass_mkf_[1] << endl << "m2jave " << d_mass_mkf_[0] << "\t" << d_mass_mkf_[1] << "\t" << d_mass_mkf_[2] << "\t" << d_mass_mkf_[3] << endl;
		out_res << endl;
//		cout << "chi_values:\t" << chi_square_q_mkf_ << "\t" << chi_square_d_mkf_ << endl;
		
	}
	outf->cd();
	outTree->Write();
	outf->Close();	
}
