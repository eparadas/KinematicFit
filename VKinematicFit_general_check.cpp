#include <iostream>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>

#include "./flat_Class.C"

void VKinematicFit_general_check()
{
	
	flat_Class t_k;
	
	
	TFile *inf = new TFile("vkinematicfit_1200_300_120.root");
	TTree *t1 = (TTree*)inf->Get("outTree");
	
	float converge[2], start_pt[8], intmd_pt[8], final_pt[8], quartet_ini_mass[2], quartet_int_mass[2], quartet_mass[2], doublet_ini_mass[4], doublet_int_mass[4], doublet_mass[4];
	int ent_mkf;
	
	TH1F h_q_ini_mass("h_q_ini_mass","h_q_ini_mass", 80, 0, 2000);
	TH1F h_q_int_mass("h_q_int_mass","h_q_int_mass", 80, 0, 2000);
	TH1F h_q_mass("h_q_mass","h_q_mass", 80, 0, 2000);
	TH1F h_q_k_mass("h_q_k_mass","h_q_k_mass", 80, 0, 2000);
	
	TH1F h_d_ini_mass("h_d_ini_mass","h_d_ini_mass", 80, 0, 2000);
	TH1F h_d_int_mass("h_d_int_mass","h_d_int_mass", 80, 0, 2000);
	TH1F h_d_mass("h_d_mass","h_d_mass", 80, 0, 2000);
	TH1F h_d_k_mass("h_d_k_mass","h_d_k_mass", 80, 0, 2000);
	
	TH2F h_pt_ini_pc_int("h_pt_ini_pc_int", "h_pt_ini_pc_int", 80, 0, 2000, 1000, -5, 5);
	TH2F h_pt_int_pc_fin("h_pt_int_pc_fin", "h_pt_int_pc_fin", 80, 0, 2000, 1000, -5, 5);
	
	t1->SetBranchAddress("event_mkf", &ent_mkf);
	t1->SetBranchAddress("converge_mkf", &converge);
	t1->SetBranchAddress("start_pt_mkf", &start_pt);
	t1->SetBranchAddress("intmd_pt_mkf", &intmd_pt);
	t1->SetBranchAddress("final_pt_mkf", &final_pt);
	t1->SetBranchAddress("quartet_ini_mass_mkf", &quartet_ini_mass);
	t1->SetBranchAddress("quartet_int_mass_mkf", &quartet_int_mass);
	t1->SetBranchAddress("quartet_mass_mkf", &quartet_mass);
	t1->SetBranchAddress("doublet_ini_mass_mkf", &doublet_ini_mass);
	t1->SetBranchAddress("doublet_int_mass_mkf", &doublet_int_mass);
	t1->SetBranchAddress("doublet_mass_mkf", &doublet_mass);
	
	for(int evt=0; evt<t1->GetEntries(); evt++)
	{
		
		t1->GetEntry(evt);
		
		t_k.GetEntry(ent_mkf);
		
		if ( converge[0] == 0 )
		{
			for(int j=0; j<2; j++)
			{
				h_q_ini_mass.Fill(quartet_ini_mass[j]);
				h_q_int_mass.Fill(quartet_int_mass[j]);
			}
			
			for(int i=0; i<2; i++)
			{
				TLorentzVector col(0,0,0,0);
				for(int j=0; j<4; j++)
				{
					TLorentzVector temp_jet(0,0,0,0);
					int temp_j = t_k.index4J[i][j];
					temp_jet.SetPtEtaPhiM( (*t_k.pt)[temp_j], (*t_k.eta)[temp_j], (*t_k.phi)[temp_j], (*t_k.mass)[temp_j] );
					col += temp_jet;
				}
				h_q_k_mass.Fill(col.M());
			}
			
			for(int j=0; j<8; j++)
				h_pt_ini_pc_int.Fill( start_pt[j], (start_pt[j] - intmd_pt[j])/start_pt[j] );

		}
		
		if ( converge[1] == 0)
		{
			h_q_mass.Fill(quartet_mass[0]);
			h_q_mass.Fill(quartet_mass[1]);
			
			for(int j=0; j<4; j++)
			{
				h_d_ini_mass.Fill(doublet_ini_mass[j]);
				h_d_int_mass.Fill(doublet_int_mass[j]);
				h_d_mass.Fill(doublet_mass[j]);
			}
			
			for(int i=0; i<4; i++)
			{
				TLorentzVector hpi(0,0,0,0);
				for(int j=0; j<2; j++)
				{
					TLorentzVector temp_jet(0,0,0,0);
					int temp_j = t_k.index2J[i][j];
					temp_jet.SetPtEtaPhiM( (*t_k.pt)[temp_j], (*t_k.eta)[temp_j], (*t_k.phi)[temp_j], (*t_k.mass)[temp_j] );
					hpi += temp_jet;
				}
				h_d_k_mass.Fill(hpi.M());
			}
			
			for(int j=0; j<8; j++)
				h_pt_int_pc_fin.Fill( intmd_pt[j], (intmd_pt[j] - final_pt[j])/intmd_pt[j] );
		}
		
	}
	
	h_q_ini_mass.SetLineColor(kBlue);
	h_q_int_mass.SetLineColor(kMagenta);
	h_q_mass.SetLineColor(kOrange);
	h_q_k_mass.SetLineColor(kGreen);
	
	h_d_ini_mass.SetLineColor(kBlue);
	h_d_int_mass.SetLineColor(kMagenta);
	h_d_mass.SetLineColor(kOrange);
	h_d_k_mass.SetLineColor(kGreen);
	
	TCanvas *c_q_mass = new TCanvas("c_q_mass", "c_q_mass", 1024, 768);
	c_q_mass->cd();
	h_q_int_mass.DrawNormalized("",1);
	h_q_ini_mass.DrawNormalized("same", 1);
	h_q_mass.DrawNormalized("same", 1);
	h_q_k_mass.DrawNormalized("same", 1);
	
	TCanvas *c_d_mass = new TCanvas("c_d_mass", "c_d_mass", 1024, 768);
	c_d_mass->cd();
	h_d_mass.DrawNormalized("",1);
	h_d_int_mass.DrawNormalized("same", 1);
	h_d_ini_mass.DrawNormalized("same", 1);
	h_d_k_mass.DrawNormalized("same", 1);
	
	TCanvas *c_pt_ini_pc_int = new TCanvas("c_pt_ini_pc_int","c_pt_ini_pc_int", 1024, 768);
	c_pt_ini_pc_int->cd();
	h_pt_ini_pc_int.Draw();
	
	TCanvas *c_pt_int_pc_fin = new TCanvas("c_pt_int_pc_fin","c_pt_int_pc_fin", 1024, 768);
	c_pt_int_pc_fin->cd();
	h_pt_int_pc_fin.Draw("colz");
	
	TFile *outf = new TFile("VKinematicFit_general_check.root", "recreate");
	h_pt_ini_pc_int.Write();
	h_pt_int_pc_fin.Write();
}
	