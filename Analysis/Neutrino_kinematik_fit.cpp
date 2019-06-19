#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TEventList.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLorentzVector.h>

using namespace std;

void Neutrino_kinematik_fit()
{
	int E_cm = 500; //CMSEnergy of samples
	int nbins = 200; //number of bins on each axis in 2D histigrams
	int cut_error = 1000; //CutError on forming equation for neutrino momenta
	string file_dir = "path/to/input/root/files/";
	string file_name = "name_of_root_file";
	TFile *file =  TFile::Open((file_dir + file_name + ".root").c_str(), "READONLY");

	std::string input_tree_name = "SLDAnalysisTree";
	string plots_dir = "/path/to/outputs/";
	TTree *tree = (TTree*)file->Get(input_tree_name.c_str());

	int r_nEntries = tree->GetEntries();
	int r_nBHadSLDecay=0;
	int r_nCHadSLDecay=0;
	int r_nSLDecayTotal=0;
	float r_EnergyCOM=0;
	double r_mcEnergyENu=0;
	double r_mcEnergyELep=0;
	double r_E_vis=0;
	double r_px_vis=0;
	double r_py_vis=0;
	double r_pz_vis=0;
	double r_E_had=0;
	double r_px_had=0;
	double r_py_had=0;
	double r_pz_had=0;
	double r_SLD_vertex_x=0;
	double r_SLD_vertex_y=0;
	double r_SLD_vertex_z=0;
	double r_Hadron_mass=0;
	double r_Mass_vis=0;
	double r_mcNu_px=0;
	double r_mcNu_py=0;
	double r_mcNu_pz=0;

	tree->SetBranchAddress("nBHadSLDecay", &r_nBHadSLDecay);
	tree->SetBranchAddress("nCHadSLDecay", &r_nCHadSLDecay);
	tree->SetBranchAddress("nSLDecayTotal", &r_nSLDecayTotal);
	tree->SetBranchAddress("EnergyCOM", &r_EnergyCOM);
	tree->SetBranchAddress("mcEnergyENu", &r_mcEnergyENu);
	tree->SetBranchAddress("mcNu_px", &r_mcNu_px);
	tree->SetBranchAddress("mcNu_py", &r_mcNu_py);
	tree->SetBranchAddress("mcNu_pz", &r_mcNu_pz);
	tree->SetBranchAddress("mcEnergyELep", &r_mcEnergyELep);
	tree->SetBranchAddress("E_vis", &r_E_vis);
	tree->SetBranchAddress("px_vis", &r_px_vis);
	tree->SetBranchAddress("py_vis", &r_py_vis);
	tree->SetBranchAddress("pz_vis", &r_pz_vis);
	tree->SetBranchAddress("Hadron_E", &r_E_had);
	tree->SetBranchAddress("Hadron_px", &r_px_had);
	tree->SetBranchAddress("Hadron_py", &r_py_had);
	tree->SetBranchAddress("Hadron_pz", &r_pz_had);
	tree->SetBranchAddress("SLD_vertex_x", &r_SLD_vertex_x);
	tree->SetBranchAddress("SLD_vertex_y", &r_SLD_vertex_y);
	tree->SetBranchAddress("SLD_vertex_z", &r_SLD_vertex_z);
	tree->SetBranchAddress("Hadron_mass", &r_Hadron_mass);
	tree->SetBranchAddress("Mass_vis", &r_Mass_vis);

	int w_nBHadSLDecay=0;
	int w_nCHadSLDecay=0;
	int w_nSLDecayTotal=0;
	float w_EnergyCOM=0;
	double w_mcEnergyENu=0;
	double w_mcEnergyELep=0;
	double w_E_vis=0;
	double w_px_vis=0;
	double w_py_vis=0;
	double w_pz_vis=0;
	double w_E_had=0;
	double w_px_had=0;
	double w_py_had=0;
	double w_pz_had=0;
	double w_SLD_vertex_x=0;
	double w_SLD_vertex_y=0;
	double w_SLD_vertex_z=0;
	double w_Hadron_mass=0;
	double w_Mass_vis=0;
	double w_mcNu_px=0;
	double w_mcNu_py=0;
	double w_mcNu_pz=0;
	int w_nEvt=0;
	double w_a_had = 0;
	double w_b_had = 0;
	double w_c_had = 0;
	double w_had_eq_err = 0;
	double w_a_nu = 0;
	double w_b_nu = 0;
	double w_c_nu = 0;
	double w_nu_eq_err = 0;
	double w_mcNeutrino_p_par = 0;
	double w_mcNeutrino_p_nor = 0;
	double w_mcNeutrino_p_par_squared = 0;
	double w_mcNeutrino_p_nor_squared = 0;
	double w_parent_p = 0;
	double w_visible_p = 0;
	double w_visible_p_par = 0;
	double w_visible_p_nor = 0;
	double w_mcNeutrino_p = 0;
	double w_mcNeutrino_p_squared = 0;
	double w_parent_p_squared = 0;
	double w_visible_p_squared = 0;
	double w_parent_p_dot_visible_p = 0;
	double w_recNeutrino_p_par_plus = 0;
	double w_recNeutrino_p_par_minus = 0;
	double w_recNeutrino_p_nor = 0;
	double w_recNeutrino_p_plus = 0;
	double w_recNeutrino_p_minus = 0;
	double w_recEnergyNu_plus = 0;
	double w_recEnergyNu_minus = 0;
	double w_recEnergyNu_close = 0;


	TFile *w_pTFile = new TFile((plots_dir+"Name_of_output_nTuple"+to_string(E_cm)+".root").c_str(), "recreate");
	TTree *w_pTTree = new TTree("NeutrinoCorrectionTree", "NeutrinoCorrectionTree");
	w_pTTree->SetDirectory(w_pTFile);
	w_pTTree->Branch("event", &w_nEvt, "event/I");
	w_pTTree->Branch("EnergyCOM", &w_EnergyCOM, "EnergyCOM/F");
	w_pTTree->Branch("nBHadSLDecay", &w_nBHadSLDecay, "nBHadSLDecay/I");
	w_pTTree->Branch("nCHadSLDecay", &w_nCHadSLDecay, "nCHadSLDecay/I");
	w_pTTree->Branch("nSLDecayTotal", &w_nSLDecayTotal, "nSLDecayTotal/I");
	w_pTTree->Branch("SLD_vertex_x", &w_SLD_vertex_x, "SLD_vertex_x/D");
	w_pTTree->Branch("SLD_vertex_y", &w_SLD_vertex_y, "SLD_vertex_y/D");
	w_pTTree->Branch("SLD_vertex_z", &w_SLD_vertex_z, "SLD_vertex_z/D");
	w_pTTree->Branch("Hadron_E", &w_E_had, "Hadron_E/D");
	w_pTTree->Branch("Hadron_px", &w_px_had, "Hadron_px/D");
	w_pTTree->Branch("Hadron_py", &w_py_had, "Hadron_py/D");
	w_pTTree->Branch("Hadron_pz", &w_pz_had, "Hadron_pz/D");
	w_pTTree->Branch("Hadron_mass", &w_Hadron_mass, "Hadron_mass/D");
	w_pTTree->Branch("E_vis", &w_E_vis, "E_vis/D");
	w_pTTree->Branch("px_vis", &w_px_vis, "px_vis/D");
	w_pTTree->Branch("py_vis", &w_py_vis, "py_vis/D");
	w_pTTree->Branch("pz_vis", &w_pz_vis, "pz_vis/D");
	w_pTTree->Branch("Mass_vis", &w_Mass_vis, "Mass_vis/D");
	w_pTTree->Branch("mcEnergyENu", &w_mcEnergyENu, "mcEnergyENu/D");
	w_pTTree->Branch("mcNu_px", &w_mcNu_px, "mcNu_px/D");
	w_pTTree->Branch("mcNu_py", &w_mcNu_py, "mcNu_py/D");
	w_pTTree->Branch("mcNu_pz", &w_mcNu_pz, "mcNu_pz/D");
	w_pTTree->Branch("mcEnergyELep", &w_mcEnergyELep, "mcEnergyELep/D");
	w_pTTree->Branch("a_had", &w_a_had, "a_had/D");
	w_pTTree->Branch("b_had", &w_b_had, "b_had/D");
	w_pTTree->Branch("c_had", &w_c_had, "c_had/D");
	w_pTTree->Branch("had_eq_err", &w_had_eq_err, "had_eq_err/D");
	w_pTTree->Branch("a_nu", &w_a_nu, "a_nu/D");
	w_pTTree->Branch("b_nu", &w_b_nu, "b_nu/D");
	w_pTTree->Branch("c_nu", &w_c_nu, "c_nu/D");
	w_pTTree->Branch("nu_eq_err", &w_nu_eq_err, "nu_eq_err/D");
	w_pTTree->Branch("mcNeutrino_p_par", &w_mcNeutrino_p_par, "mcNeutrino_p_par/D");
	w_pTTree->Branch("mcNeutrino_p_nor", &w_mcNeutrino_p_nor, "mcNeutrino_p_nor/D");
	w_pTTree->Branch("mcNeutrino_p_par_squared", &w_mcNeutrino_p_par_squared, "mcNeutrino_p_par_squared/D");
	w_pTTree->Branch("mcNeutrino_p_nor_squared", &w_mcNeutrino_p_nor_squared, "mcNeutrino_p_nor_squared/D");
	w_pTTree->Branch("parent_p", &w_parent_p, "parent_p/D");
	w_pTTree->Branch("visible_p", &w_visible_p, "visible_p/D");
	w_pTTree->Branch("visible_p_par", &w_visible_p_par, "visible_p_par/D");
	w_pTTree->Branch("visible_p_nor", &w_visible_p_nor, "visible_p_nor/D");
	w_pTTree->Branch("mcNeutrino_p", &w_mcNeutrino_p, "mcNeutrino_p/D");
	w_pTTree->Branch("mcNeutrino_p_squared", &w_mcNeutrino_p_squared, "mcNeutrino_p_squared/D");
	w_pTTree->Branch("parent_p_squared", &w_parent_p_squared, "parent_p_squared/D");
	w_pTTree->Branch("visible_p_squared", &w_visible_p_squared, "visible_p_squared/D");
	w_pTTree->Branch("parent_p_dot_visible_p", &w_parent_p_dot_visible_p, "parent_p_dot_visible_p/D");
	w_pTTree->Branch("recNeutrino_p_par_plus", &w_recNeutrino_p_par_plus, "recNeutrino_p_par_plus/D");
	w_pTTree->Branch("recNeutrino_p_par_minus", &w_recNeutrino_p_par_minus, "recNeutrino_p_par_minus/D");
	w_pTTree->Branch("recNeutrino_p_nor", &w_recNeutrino_p_nor, "recNeutrino_p_nor/D");
	w_pTTree->Branch("recNeutrino_p_plus", &w_recNeutrino_p_plus, "recNeutrino_p_plus/D");
	w_pTTree->Branch("recNeutrino_p_minus", &w_recNeutrino_p_minus, "recNeutrino_p_minus/D");
	w_pTTree->Branch("recEnergyNu_plus", &w_recEnergyNu_plus, "recEnergyNu_plus/D");
	w_pTTree->Branch("recEnergyNu_minus", &w_recEnergyNu_minus, "recEnergyNu_minus/D");
	w_pTTree->Branch("recEnergyNu_close", &w_recEnergyNu_close, "recEnergyNu_close/D");

	TH2D *Pnu_PbPvis = new TH2D("Pnu_PbPvis","; p_{#nu,#parallel} [GeV]; p_{Had} - p_{vis,#parallel} [GeV]", nbins, 0., 200, nbins, 0., 200 );
	Pnu_PbPvis->GetXaxis()->SetTitleSize(0.05); Pnu_PbPvis->GetYaxis()->SetTitleSize(0.05); Pnu_PbPvis->GetXaxis()->SetLabelSize(0.04); Pnu_PbPvis->GetYaxis()->SetLabelSize(0.04);

	TH2D *Pnu_Pvis = new TH2D("Pnu_Pvis","; p_{#nu,#perp} [GeV]; p_{vis,#perp} [GeV]", nbins, 0., 4., nbins, 0., 4. );
	Pnu_Pvis->GetXaxis()->SetTitleSize(0.05); Pnu_Pvis->GetYaxis()->SetTitleSize(0.05); Pnu_Pvis->GetXaxis()->SetLabelSize(0.04); Pnu_Pvis->GetYaxis()->SetLabelSize(0.04);

	TH2D *Enu_EbEvis = new TH2D("Enu_EbEvis","; E_{#nu} [GeV]; E_{Had} - E_{vis} [GeV]", nbins, 0., 200, nbins, 0., 200 );
	Enu_EbEvis->GetXaxis()->SetTitleSize(0.05); Enu_EbEvis->GetYaxis()->SetTitleSize(0.05); Enu_EbEvis->GetXaxis()->SetLabelSize(0.04); Enu_EbEvis->GetYaxis()->SetLabelSize(0.04);

	TH2D *Pnu_Enu = new TH2D("Pnu_Enu","; p_{#nu} [GeV]; E_{#nu} [GeV]", nbins, 0., 200, nbins, 0., 200 );
	Pnu_Enu->GetXaxis()->SetTitleSize(0.05); Pnu_Enu->GetYaxis()->SetTitleSize(0.05); Pnu_Enu->GetXaxis()->SetLabelSize(0.04); Pnu_Enu->GetYaxis()->SetLabelSize(0.04);

	TH2D *Pnu2_PbPvis2 = new TH2D("Pnu2_PbPvis2","; p_{#nu}^{2} [GeV^{2}]; p_{Had}^{2} + p_{vis}^{2} - 2#vec{p}_{Had}.#vec{p}_{vis} [GeV^{2}]", nbins, 0., 40000, nbins, 0., 40000 );
	Pnu2_PbPvis2->GetXaxis()->SetTitleSize(0.06); Pnu2_PbPvis2->GetYaxis()->SetTitleSize(0.06); Pnu2_PbPvis2->GetXaxis()->SetLabelSize(0.03); Pnu2_PbPvis2->GetYaxis()->SetLabelSize(0.03);

	TH2D *Enu2_EbEvis2 = new TH2D("Enu2_EbEvis2","; E_{#nu}^{2} [GeV^{2}]; (E_{Had} - E_{vis})^{2} [GeV^{2}]", nbins, 0., 40000, nbins, 0., 40000 );
	Enu2_EbEvis2->GetXaxis()->SetTitleSize(0.06); Enu2_EbEvis2->GetYaxis()->SetTitleSize(0.06); Enu2_EbEvis2->GetXaxis()->SetLabelSize(0.03); Enu2_EbEvis2->GetYaxis()->SetLabelSize(0.03);

	TH2D *EbEvis2_PbPvis2 = new TH2D("EbEvis2_PbPvis2","; E_{Had}^{2} + E_{vis}^{2} - 2 E_{Had}.E_{vis} [GeV^{2}]; p_{Had}^{2} + p_{vis}^{2} - 2 p_{Had}.p_{vis,#parallel} [GeV^{2}]", nbins, 0., 40000, nbins, 0., 40000 );
	EbEvis2_PbPvis2->GetXaxis()->SetTitleSize(0.06); EbEvis2_PbPvis2->GetYaxis()->SetTitleSize(0.06); EbEvis2_PbPvis2->GetXaxis()->SetLabelSize(0.03); EbEvis2_PbPvis2->GetYaxis()->SetLabelSize(0.03);

	TH2D *mbmvis2_PbPvis2 = new TH2D("mbmvis2_PbPvis2","; m_{Had}^{2} + m_{vis}^{2} + 2 p_{Had}.p_{vis,#parallel} [GeV^{2}]; 2 E_{vis}.#sqrt{m_{Had}^{2} + p_{Had}^{2}} [GeV^{2}]", nbins, 0., 140000, nbins, 0., 140000 );
	mbmvis2_PbPvis2->GetXaxis()->SetTitleSize(0.06); mbmvis2_PbPvis2->GetYaxis()->SetTitleSize(0.06); mbmvis2_PbPvis2->GetXaxis()->SetLabelSize(0.03); mbmvis2_PbPvis2->GetYaxis()->SetLabelSize(0.03);

	TH2D *mbmvis2_PbPvis2_squared = new TH2D("mbmvis2_PbPvis2_squared","; (m_{Had}^{2} + m_{vis}^{2})^{2} + 4 p_{Had}^{2}.p_{vis,#parallel}^{2} + 4 p_{Had}.p_{vis,#parallel}(m_{Had}^{2} + m_{vis}^{2}) [GeV^{4}]; 4 E_{vis}^{2} (m_{Had}^{2} + p_{Had}^{2}) [GeV^{4}]", nbins, 0., 20000000000, nbins, 0., 20000000000 );
	mbmvis2_PbPvis2_squared->GetXaxis()->SetTitleSize(0.045); mbmvis2_PbPvis2_squared->GetYaxis()->SetTitleSize(0.05); mbmvis2_PbPvis2_squared->GetXaxis()->SetLabelSize(0.03); mbmvis2_PbPvis2_squared->GetYaxis()->SetLabelSize(0.03);

	TH1D* Pb_eqerr = new TH1D("Pb_eqerr","4(p_{vis,#parallel}^{2}-E_{vis}^{2})p_{Had}^{2} + 4p_{vis,#parallel}(m_{Had}^{2}+m_{vis}^{2})p_{Had} + (m_{Had}^{2}+m_{vis}^{2})^{2} - 4m_{Had}^{2}E_{vis}^{2} [GeV^{4}]", 2000,-1000,1000);

	TH1D* Pnu_eqerr = new TH1D("Pnu_eqerr","4(p_{vis,#parallel}^{2}-E_{vis}^{2})p_{#nu,#parallel}^{2} + 4p_{vis,#parallel}(m_{Had}^{2} - m_{vis}^{2} - 2p_{vis,#perp}^{2})p_{#nu,#parallel} + (m_{Had}^{2}+m_{vis}^{2})^{2} - 4m_{Had}^{2}E_{vis}^{2} + 4p_{vis,#parallel}(m_{Had}^{2} - p_{vis,#perp}^{2})", 2000,-1000.,1000.);

	TH2D *mcENu_recENu_plus = new TH2D("mcENu_recENu_plus","; E_{#nu}^{MC} [GeV]; E_{#nu}^{REC+} [GeV]", nbins, 0., 200., nbins, 0., 200. );
	mcENu_recENu_plus->GetXaxis()->SetTitleSize(0.05); mcENu_recENu_plus->GetYaxis()->SetTitleSize(0.05); mcENu_recENu_plus->GetXaxis()->SetLabelSize(0.04); mcENu_recENu_plus->GetYaxis()->SetLabelSize(0.04);

	TH2D *mcENu_recENu_minus = new TH2D("mcENu_recENu_minus","; E_{#nu}^{MC} [GeV]; E_{#nu}^{REC-} [GeV]", nbins, 0., 200., nbins, 0., 200. );
	mcENu_recENu_minus->GetXaxis()->SetTitleSize(0.05); mcENu_recENu_minus->GetYaxis()->SetTitleSize(0.05); mcENu_recENu_minus->GetXaxis()->SetLabelSize(0.04); mcENu_recENu_minus->GetYaxis()->SetLabelSize(0.04);

	TH2D *mcENu_recENu_close = new TH2D("mcENu_recENu_close","; E_{#nu}^{MC} [GeV]; E_{#nu}^{REC+/-} [GeV]", nbins, 0., 200., nbins, 0., 200. );
	mcENu_recENu_close->GetXaxis()->SetTitleSize(0.05); mcENu_recENu_close->GetYaxis()->SetTitleSize(0.05); mcENu_recENu_close->GetXaxis()->SetLabelSize(0.04); mcENu_recENu_close->GetYaxis()->SetLabelSize(0.04);

	for (int i = 0; i<r_nEntries; ++i)
	{
		tree->GetEntry(i);
		w_nEvt = i;
		w_nBHadSLDecay = r_nBHadSLDecay;
		w_nCHadSLDecay = r_nCHadSLDecay;
		w_nSLDecayTotal = r_nBHadSLDecay + r_nCHadSLDecay;
		w_SLD_vertex_x = r_SLD_vertex_x;
		w_SLD_vertex_y = r_SLD_vertex_y;
		w_SLD_vertex_z = r_SLD_vertex_z;
		w_EnergyCOM = r_EnergyCOM;

		TVector3 parent_p(r_px_had,r_py_had,r_pz_had);
		TLorentzVector parent_tlv(parent_p,r_E_had);
		double parent_mass = parent_tlv.M();
		w_E_had = parent_tlv.E();
		w_px_had = parent_tlv.Px();
		w_py_had = parent_tlv.Py();
		w_pz_had = parent_tlv.Pz();
		w_Hadron_mass = parent_mass;
		w_parent_p = parent_p.Mag();
		w_parent_p_squared = parent_p.Mag2();

		TVector3 visible_p(r_px_vis,r_py_vis,r_pz_vis);
		TLorentzVector visible_tlv(visible_p,r_E_vis);
		double visible_mass = visible_tlv.M();
		w_E_vis = visible_tlv.E();
		w_px_vis = visible_tlv.Px();
		w_py_vis = visible_tlv.Py();
		w_pz_vis = visible_tlv.Pz();
		w_Mass_vis = visible_mass;
		w_visible_p = visible_p.Mag();
		w_visible_p_squared = visible_p.Mag2();


		TVector3 mcNeutrino_p(r_mcNu_px,r_mcNu_py,r_mcNu_pz);
		TLorentzVector mcNeutrino_tlv(mcNeutrino_p,r_mcEnergyENu);
		w_mcEnergyENu = r_mcEnergyENu;
		w_mcNu_px = r_mcNu_px;
		w_mcNu_py = r_mcNu_py;
		w_mcNu_pz = r_mcNu_pz;

		TVector3 u(parent_p.Unit());
		TVector3 parent_p_par(parent_p.Dot(u)*u);
		TVector3 parent_p_nor(parent_p - parent_p_par);
		TVector3 visible_p_par(visible_p.Dot(u)*u);
		TVector3 visible_p_nor(visible_p - visible_p_par);
		TVector3 mcNeutrino_p_par(mcNeutrino_p.Dot(u)*u);
		TVector3 mcNeutrino_p_nor(mcNeutrino_p - mcNeutrino_p_par);

		double a_had = 4*(visible_p_par.Mag2() - pow(visible_tlv.E(),2));
		double b_had = 4*visible_p_par.Mag()*(parent_tlv.M2() + visible_tlv.M2());
		double c_had = pow(parent_tlv.M2() + visible_tlv.M2(),2) - 4*parent_tlv.M2()*pow(visible_tlv.E(),2);
		double had_eq_err = a_had * parent_p.Mag2() + b_had * parent_p.Mag() + c_had;
		w_a_had = a_had;
		w_b_had = b_had;
		w_c_had = c_had;
		w_had_eq_err = had_eq_err;

		double a_nu = 4*(visible_p_par.Mag2() - pow(visible_tlv.E(),2));
		double b_nu = 4*visible_p_par.Mag()*(parent_tlv.M2() - visible_tlv.M2() - 2*visible_p_nor.Mag2());
		double c_nu = pow(parent_tlv.M2() - visible_tlv.M2(),2) - 4*visible_p_nor.Mag2()*(parent_tlv.M2() + visible_p_par.Mag2());
		double nu_eq_err = a_nu * mcNeutrino_p.Mag2() + b_nu * mcNeutrino_p.Mag() + c_nu; // "a_nu * mcNeutrino_p.Mag2() + b_nu * mcNeutrino_p.Mag() + c_nu" should be 0 in second hand but is not!
		w_a_nu = a_nu;
		w_b_nu = b_nu;
		w_c_nu = c_nu;
		w_nu_eq_err = nu_eq_err;

		w_visible_p_par = visible_p_par.Mag();
		w_visible_p_nor = visible_p_nor.Mag();

		w_parent_p_dot_visible_p = parent_p.Dot(visible_p);

		w_mcNeutrino_p_par = mcNeutrino_p_par.Mag();
		w_mcNeutrino_p_nor = mcNeutrino_p_nor.Mag();
		w_mcNeutrino_p_par_squared = mcNeutrino_p_par.Mag2();
		w_mcNeutrino_p_nor_squared = mcNeutrino_p_nor.Mag2();
		w_mcNeutrino_p = mcNeutrino_p.Mag();
		w_mcNeutrino_p_squared = mcNeutrino_p.Mag2();

		TVector3 recNeutrino_p_par_plus((-b_nu - sqrt(pow(b_nu,2) - 4 * a_nu * c_nu))/(2 * a_nu) * u);
		TVector3 recNeutrino_p_par_minus((-b_nu + sqrt(pow(b_nu,2) - 4 * a_nu * c_nu))/(2 * a_nu) * u);
		TVector3 recNeutrino_p_nor(-1*visible_p_nor);
		TVector3 recNeutrino_p_plus(recNeutrino_p_par_plus + recNeutrino_p_nor);
		TVector3 recNeutrino_p_minus(recNeutrino_p_par_minus + recNeutrino_p_nor);

		w_recNeutrino_p_par_plus = recNeutrino_p_par_plus.Mag();
		w_recNeutrino_p_par_minus = recNeutrino_p_par_minus.Mag();
		w_recNeutrino_p_nor = recNeutrino_p_nor.Mag();
		w_recNeutrino_p_plus = recNeutrino_p_plus.Mag();
		w_recNeutrino_p_minus = recNeutrino_p_minus.Mag();


		double recEnergyNu_plus = recNeutrino_p_plus.Mag();
		double recEnergyNu_minus = recNeutrino_p_minus.Mag();
		double recEnergyNu_close = 0;
		if (abs(recEnergyNu_plus - r_mcEnergyENu) < abs(recEnergyNu_minus - r_mcEnergyENu))
		{
			recEnergyNu_close = recEnergyNu_plus;
		}
		else
		{
			recEnergyNu_close = recEnergyNu_minus;
		}
		w_recEnergyNu_plus = recEnergyNu_plus;
		w_recEnergyNu_minus = recEnergyNu_minus;
		w_recEnergyNu_close = recEnergyNu_close;

		w_pTTree->Fill();

		if (abs(nu_eq_err) >= cut_error) // fill histograms if the error is less than CutError
			continue;

		if (r_nBHadSLDecay + r_nCHadSLDecay == 1) // Perform (plot) netrino correction for events with only ONE semi-leptonic decay
		{

			Pnu_PbPvis->Fill(mcNeutrino_p_par.Mag(),(parent_p - visible_p_par).Mag());
			Pnu_Pvis->Fill(mcNeutrino_p_nor.Mag(),visible_p_nor.Mag());
			Enu_EbEvis->Fill(mcNeutrino_tlv.E(),parent_tlv.E() - visible_tlv.E());
			Pnu_Enu->Fill(mcNeutrino_p.Mag(),mcNeutrino_tlv.E());
			Pnu2_PbPvis2->Fill(mcNeutrino_p.Mag2(),parent_p.Mag2()+visible_p.Mag2()-2*parent_p.Dot(visible_p));
			Enu2_EbEvis2->Fill(pow(mcNeutrino_tlv.E(),2),pow(parent_tlv.E()-visible_tlv.E(),2));
			EbEvis2_PbPvis2->Fill(pow(parent_tlv.E(),2)+pow(visible_tlv.E(),2)-2*parent_tlv.E()*visible_tlv.E(),parent_p.Mag2()+visible_p.Mag2()-2*parent_p.Mag()*visible_p_par.Mag());
			mbmvis2_PbPvis2->Fill(parent_tlv.M2()+visible_tlv.M2()+2*parent_p.Mag()*visible_p_par.Mag(),2*visible_tlv.E()*sqrt(parent_tlv.M2()+parent_p.Mag2()));
			mbmvis2_PbPvis2_squared->Fill(pow(parent_tlv.M2()+visible_tlv.M2()+2*parent_p.Mag()*visible_p_par.Mag(),2),pow(2*visible_tlv.E()*sqrt(parent_tlv.M2()+parent_p.Mag2()),2));
			Pb_eqerr->Fill(had_eq_err);
			Pnu_eqerr->Fill(nu_eq_err);
			mcENu_recENu_plus->Fill(r_mcEnergyENu,recEnergyNu_plus);
			mcENu_recENu_minus->Fill(r_mcEnergyENu,recEnergyNu_minus);
			mcENu_recENu_close->Fill(r_mcEnergyENu,recEnergyNu_close);
		}
	}

	// plot histograms

	TCanvas *can_mcENu_recENu_plus = new TCanvas("can_mcENu_recENu_plus", "mcENu_recENu_plus", 1280,1000);
	mcENu_recENu_plus->Draw("colz");
	gPad->Update();
	can_mcENu_recENu_plus->Modified();
	can_mcENu_recENu_plus->Update();
	TPaveStats *tps_mcENu_recENu_plus = (TPaveStats*)mcENu_recENu_plus->FindObject("stats");
	TPaletteAxis *pal_mcENu_recENu_plus = (TPaletteAxis*)mcENu_recENu_plus->GetListOfFunctions()->FindObject("palette");
	can_mcENu_recENu_plus->Update();
	tps_mcENu_recENu_plus->SetX1NDC(0.65);
	tps_mcENu_recENu_plus->SetY1NDC(0.2);
	tps_mcENu_recENu_plus->SetX2NDC(0.80);
	tps_mcENu_recENu_plus->SetY2NDC(0.465);
	pal_mcENu_recENu_plus->SetX1NDC(0.85);
	pal_mcENu_recENu_plus->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_mcENu_recENu_plus->Modified();
	can_mcENu_recENu_plus->Update();
	can_mcENu_recENu_plus->SaveAs((plots_dir+"mcENu_recENu_plus_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_mcENu_recENu_minus = new TCanvas("can_mcENu_recENu_minus", "mcENu_recENu_minus", 1280,1000);
	mcENu_recENu_minus->Draw("colz");
	gPad->Update();
	can_mcENu_recENu_minus->Modified();
	can_mcENu_recENu_minus->Update();
	TPaveStats *tps_mcENu_recENu_minus = (TPaveStats*)mcENu_recENu_minus->FindObject("stats");
	TPaletteAxis *pal_mcENu_recENu_minus = (TPaletteAxis*)mcENu_recENu_minus->GetListOfFunctions()->FindObject("palette");
	can_mcENu_recENu_minus->Update();
	tps_mcENu_recENu_minus->SetX1NDC(0.65);
	tps_mcENu_recENu_minus->SetY1NDC(0.2);
	tps_mcENu_recENu_minus->SetX2NDC(0.80);
	tps_mcENu_recENu_minus->SetY2NDC(0.465);
	pal_mcENu_recENu_minus->SetX1NDC(0.85);
	pal_mcENu_recENu_minus->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_mcENu_recENu_minus->Modified();
	can_mcENu_recENu_minus->Update();
	can_mcENu_recENu_minus->SaveAs((plots_dir+"mcENu_recENu_minus_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_mcENu_recENu_close = new TCanvas("can_mcENu_recENu_close", "mcENu_recENu_close", 1280,1000);
	mcENu_recENu_close->Draw("colz");
	gPad->Update();
	can_mcENu_recENu_close->Modified();
	can_mcENu_recENu_close->Update();
	TPaveStats *tps_mcENu_recENu_close = (TPaveStats*)mcENu_recENu_close->FindObject("stats");
	TPaletteAxis *pal_mcENu_recENu_close = (TPaletteAxis*)mcENu_recENu_close->GetListOfFunctions()->FindObject("palette");
	can_mcENu_recENu_close->Update();
	tps_mcENu_recENu_close->SetX1NDC(0.65);
	tps_mcENu_recENu_close->SetY1NDC(0.2);
	tps_mcENu_recENu_close->SetX2NDC(0.80);
	tps_mcENu_recENu_close->SetY2NDC(0.465);
	pal_mcENu_recENu_close->SetX1NDC(0.85);
	pal_mcENu_recENu_close->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_mcENu_recENu_close->Modified();
	can_mcENu_recENu_close->Update();
	can_mcENu_recENu_close->SaveAs((plots_dir+"mcENu_recENu_close_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());
	delete can_mcENu_recENu_close;


	TCanvas *can_Pnu_PbPvis = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Pnu_PbPvis->Draw("colz");
	gPad->Update();
	can_Pnu_PbPvis->Modified();
	can_Pnu_PbPvis->Update();
	TPaveStats *tps_Pnu_PbPvis = (TPaveStats*)Pnu_PbPvis->FindObject("stats");
	TPaletteAxis *pal_Pnu_PbPvis = (TPaletteAxis*)Pnu_PbPvis->GetListOfFunctions()->FindObject("palette");
	can_Pnu_PbPvis->Update();
	tps_Pnu_PbPvis->SetX1NDC(0.65);
	tps_Pnu_PbPvis->SetY1NDC(0.2);
	tps_Pnu_PbPvis->SetX2NDC(0.80);
	tps_Pnu_PbPvis->SetY2NDC(0.465);
	pal_Pnu_PbPvis->SetX1NDC(0.85);
	pal_Pnu_PbPvis->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Pnu_PbPvis->Modified();
	can_Pnu_PbPvis->Update();
	can_Pnu_PbPvis->SaveAs((plots_dir+"Pnu_PbPvis_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());


	TCanvas *can_Pnu_Pvis = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Pnu_Pvis->Draw("colz");
	gPad->Update();
	can_Pnu_Pvis->Modified();
	can_Pnu_Pvis->Update();
	TPaveStats *tps_Pnu_Pvis = (TPaveStats*)Pnu_Pvis->FindObject("stats");
	TPaletteAxis *pal_Pnu_Pvis = (TPaletteAxis*)Pnu_Pvis->GetListOfFunctions()->FindObject("palette");
	can_Pnu_Pvis->Update();
	tps_Pnu_Pvis->SetX1NDC(0.65);
	tps_Pnu_Pvis->SetY1NDC(0.2);
	tps_Pnu_Pvis->SetX2NDC(0.80);
	tps_Pnu_Pvis->SetY2NDC(0.465);
	pal_Pnu_Pvis->SetX1NDC(0.85);
	pal_Pnu_Pvis->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Pnu_Pvis->Modified();
	can_Pnu_Pvis->Update();
	can_Pnu_Pvis->SaveAs((plots_dir+"Pnu_Pvis_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_Enu_EbEvis = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Enu_EbEvis->Draw("colz");
	gPad->Update();
	can_Enu_EbEvis->Modified();
	can_Enu_EbEvis->Update();
	TPaveStats *tps_Enu_EbEvis = (TPaveStats*)Enu_EbEvis->FindObject("stats");
	TPaletteAxis *pal_Enu_EbEvis = (TPaletteAxis*)Enu_EbEvis->GetListOfFunctions()->FindObject("palette");
	can_Enu_EbEvis->Update();
	tps_Enu_EbEvis->SetX1NDC(0.65);
	tps_Enu_EbEvis->SetY1NDC(0.2);
	tps_Enu_EbEvis->SetX2NDC(0.80);
	tps_Enu_EbEvis->SetY2NDC(0.465);
	pal_Enu_EbEvis->SetX1NDC(0.85);
	pal_Enu_EbEvis->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Enu_EbEvis->Modified();
	can_Enu_EbEvis->Update();
	can_Enu_EbEvis->SaveAs((plots_dir+"Enu_EbEvis_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_Pnu_Enu = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Pnu_Enu->Draw("colz");
	gPad->Update();
	can_Pnu_Enu->Modified();
	can_Pnu_Enu->Update();
	TPaveStats *tps_Pnu_Enu = (TPaveStats*)Pnu_Enu->FindObject("stats");
	TPaletteAxis *pal_Pnu_Enu = (TPaletteAxis*)Pnu_Enu->GetListOfFunctions()->FindObject("palette");
	can_Pnu_Enu->Update();
	tps_Pnu_Enu->SetX1NDC(0.65);
	tps_Pnu_Enu->SetY1NDC(0.2);
	tps_Pnu_Enu->SetX2NDC(0.80);
	tps_Pnu_Enu->SetY2NDC(0.465);
	pal_Pnu_Enu->SetX1NDC(0.85);
	pal_Pnu_Enu->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Pnu_Enu->Modified();
	can_Pnu_Enu->Update();
	can_Pnu_Enu->SaveAs((plots_dir+"Pnu_Enu_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_Enu2_EbEvis2 = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Enu2_EbEvis2->Draw("colz");
	gPad->Update();
	can_Enu2_EbEvis2->Modified();
	can_Enu2_EbEvis2->Update();
	TPaveStats *tps_Enu2_EbEvis2 = (TPaveStats*)Enu2_EbEvis2->FindObject("stats");
	TPaletteAxis *pal_Enu2_EbEvis2 = (TPaletteAxis*)Enu2_EbEvis2->GetListOfFunctions()->FindObject("palette");
	can_Enu2_EbEvis2->Update();
	tps_Enu2_EbEvis2->SetX1NDC(0.65);
	tps_Enu2_EbEvis2->SetY1NDC(0.2);
	tps_Enu2_EbEvis2->SetX2NDC(0.80);
	tps_Enu2_EbEvis2->SetY2NDC(0.465);
	pal_Enu2_EbEvis2->SetX1NDC(0.85);
	pal_Enu2_EbEvis2->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Enu2_EbEvis2->Modified();
	can_Enu2_EbEvis2->Update();
	can_Enu2_EbEvis2->SaveAs((plots_dir+"Enu2_EbEvis2_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_EbEvis2_PbPvis2 = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	EbEvis2_PbPvis2->Draw("colz");
	gPad->Update();
	can_EbEvis2_PbPvis2->Modified();
	can_EbEvis2_PbPvis2->Update();
	TPaveStats *tps_EbEvis2_PbPvis2 = (TPaveStats*)EbEvis2_PbPvis2->FindObject("stats");
	TPaletteAxis *pal_EbEvis2_PbPvis2 = (TPaletteAxis*)EbEvis2_PbPvis2->GetListOfFunctions()->FindObject("palette");
	can_EbEvis2_PbPvis2->Update();
	tps_EbEvis2_PbPvis2->SetX1NDC(0.65);
	tps_EbEvis2_PbPvis2->SetY1NDC(0.2);
	tps_EbEvis2_PbPvis2->SetX2NDC(0.80);
	tps_EbEvis2_PbPvis2->SetY2NDC(0.465);
	pal_EbEvis2_PbPvis2->SetX1NDC(0.85);
	pal_EbEvis2_PbPvis2->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_EbEvis2_PbPvis2->Modified();
	can_EbEvis2_PbPvis2->Update();
	can_EbEvis2_PbPvis2->SaveAs((plots_dir+"EbEvis2_PbPvis2_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_mbmvis2_PbPvis2 = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	mbmvis2_PbPvis2->Draw("colz");
	gPad->Update();
	can_mbmvis2_PbPvis2->Modified();
	can_mbmvis2_PbPvis2->Update();
	TPaveStats *tps_mbmvis2_PbPvis2 = (TPaveStats*)mbmvis2_PbPvis2->FindObject("stats");
	TPaletteAxis *pal_mbmvis2_PbPvis2 = (TPaletteAxis*)mbmvis2_PbPvis2->GetListOfFunctions()->FindObject("palette");
	can_mbmvis2_PbPvis2->Update();
	tps_mbmvis2_PbPvis2->SetX1NDC(0.65);
	tps_mbmvis2_PbPvis2->SetY1NDC(0.2);
	tps_mbmvis2_PbPvis2->SetX2NDC(0.80);
	tps_mbmvis2_PbPvis2->SetY2NDC(0.465);
	pal_mbmvis2_PbPvis2->SetX1NDC(0.85);
	pal_mbmvis2_PbPvis2->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_mbmvis2_PbPvis2->Modified();
	can_mbmvis2_PbPvis2->Update();
	can_mbmvis2_PbPvis2->SaveAs((plots_dir+"mbmvis2_PbPvis2_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_mbmvis2_PbPvis2_squared = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	mbmvis2_PbPvis2_squared->Draw("colz");
	gPad->Update();
	can_mbmvis2_PbPvis2_squared->Modified();
	can_mbmvis2_PbPvis2_squared->Update();
	TPaveStats *tps_mbmvis2_PbPvis2_squared = (TPaveStats*)mbmvis2_PbPvis2_squared->FindObject("stats");
	TPaletteAxis *pal_mbmvis2_PbPvis2_squared = (TPaletteAxis*)mbmvis2_PbPvis2_squared->GetListOfFunctions()->FindObject("palette");
	can_mbmvis2_PbPvis2_squared->Update();
	tps_mbmvis2_PbPvis2_squared->SetX1NDC(0.65);
	tps_mbmvis2_PbPvis2_squared->SetY1NDC(0.2);
	tps_mbmvis2_PbPvis2_squared->SetX2NDC(0.80);
	tps_mbmvis2_PbPvis2_squared->SetY2NDC(0.465);
	pal_mbmvis2_PbPvis2_squared->SetX1NDC(0.85);
	pal_mbmvis2_PbPvis2_squared->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_mbmvis2_PbPvis2_squared->Modified();
	can_mbmvis2_PbPvis2_squared->Update();
	can_mbmvis2_PbPvis2_squared->SaveAs((plots_dir+"mbmvis2_PbPvis2_squared_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	TCanvas *can_Pnu2_PbPvis2 = new TCanvas("can", "mcENu_recENu_close", 1280,1000);
	Pnu2_PbPvis2->Draw("colz");
	gPad->Update();
	can_Pnu2_PbPvis2->Modified();
	can_Pnu2_PbPvis2->Update();
	TPaveStats *tps_Pnu2_PbPvis2 = (TPaveStats*)Pnu2_PbPvis2->FindObject("stats");
	TPaletteAxis *pal_Pnu2_PbPvis2 = (TPaletteAxis*)Pnu2_PbPvis2->GetListOfFunctions()->FindObject("palette");
	can_Pnu2_PbPvis2->Update();
	tps_Pnu2_PbPvis2->SetX1NDC(0.65);
	tps_Pnu2_PbPvis2->SetY1NDC(0.2);
	tps_Pnu2_PbPvis2->SetX2NDC(0.80);
	tps_Pnu2_PbPvis2->SetY2NDC(0.465);
	pal_Pnu2_PbPvis2->SetX1NDC(0.85);
	pal_Pnu2_PbPvis2->SetX2NDC(0.885);
	gPad->SetRightMargin(0.15);
	gPad->Update();
	can_Pnu2_PbPvis2->Modified();
	can_Pnu2_PbPvis2->Update();
	can_Pnu2_PbPvis2->SaveAs((plots_dir+"Pnu2_PbPvis2_bb_"+ to_string(E_cm) + "_GeV_EqErr_" + to_string(cut_error) + ".pdf").c_str());

	// save output root file

	w_pTFile->cd();
	w_pTTree->Write();
	w_pTFile->Close();
	delete w_pTFile;


}
