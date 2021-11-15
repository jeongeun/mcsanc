#include<Riostream.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"

Double_t getEKfactorMid(double x, bool IsAdditive){
	if(x > 7000.0) { x = 7000 ; }
	//Median of Add Fac Fit
	if(IsAdditive == 0) return ( 1.21213 -(0.00019292*x) + (7.02667e-08*pow(x,2)) -(1.09576e-11*pow(x,3)) + (5.17257e-16*pow(x,4)) );
	//Additive Fit
	if(IsAdditive == 1) return (1.24039  -(0.000235493*x) + (1.03222e-07*pow(x,2)) -(1.93752e-11 *pow(x,3)) + (1.18182e-15*pow(x,4)) );
}

Double_t getMKfactorMid(double x, bool IsAdditive){
	if(x > 7000.0) { x = 7000 ; }
	//Middle of Add Fac Fit
	if(IsAdditive == 0) return (  1.21384 -(0.00020136*x) + (8.32023e-08*pow(x,2)) -(1.55027e-11*pow(x,3)) + (9.41259e-16*pow(x,4)) );
	//Additive Fit
	if(IsAdditive == 1) return (1.24197 -(0.000243528*x) + (1.15169e-07*pow(x,2)) -(2.35027e-11*pow(x,3)) + (1.56328e-15*pow(x,4)) );
}

void histMake(){

   TFile* BaseLO   = new TFile("Wmass.root" , "READ");       //from Sebastian
   TH1D* h_base_lo = (TH1D*)BaseLO-> Get("Wgenmass");  

   TString datap_ewk_lo   = "./backup/LO_EW_Wp_e_m34_hist_new.dat" ;  //MCSANC LO W+ 
   TString datam_ewk_lo   = "./backup/LO_EW_Wm_e_m34_hist_new.dat" ;  //MCSANC LO W-
   TString datap_ewk_nlo  = "./backup/NLO_EW_Wp_e_m34_hist_new.dat";  //MCSANC NLO W+ 
   TString datam_ewk_nlo  = "./backup/NLO_EW_Wm_e_m34_hist_new.dat";  //MCSANC NLO W- 
   TString datap_qcd_lo   = "./backup/LO_QCD_Wp_hist.dat"          ;  //FEWZ LO W+         
   TString datam_qcd_lo   = "./backup/LO_QCD_Wm_hist.dat"          ;  //FEWZ LO W-         
   TString datap_qcd_nnlo = "./backup/NNLO_QCD_Wp_hist.dat"        ;  //FEWZ NNLO W+       
   TString datam_qcd_nnlo = "./backup/NNLO_QCD_Wm_hist.dat"        ;  //FEWZ NNLO W-       

   ifstream  inp_ewk_lo, inm_ewk_lo, inp_ewk_nlo, inm_ewk_nlo, inp_qcd_lo, inm_qcd_lo, inp_qcd_nnlo, inm_qcd_nnlo; 

   //Open Data (EWK MCSANC)
   inp_ewk_lo.open(datap_ewk_lo);
   inm_ewk_lo.open(datam_ewk_lo);
   inp_ewk_nlo.open(datap_ewk_nlo);
   inm_ewk_nlo.open(datam_ewk_nlo);

   ///Open Data (QCD FEWZ)
   inp_qcd_lo.open(datap_qcd_lo);
   inm_qcd_lo.open(datam_qcd_lo);
   inp_qcd_nnlo.open(datap_qcd_nnlo);
   inm_qcd_nnlo.open(datam_qcd_nnlo);

   TH1D *h_e_EW_LO = new TH1D("h_e_EW_LO" ,"h_e_EW_LO MCSANC ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_EW_LO ->Sumw2(); 
   TH1D *h_e_EW_NLO= new TH1D("h_e_EW_NLO","h_e_EW_NLO MCSANC; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_EW_NLO->Sumw2();
 
   TH1D *h_QCD_LO   = new TH1D("h_QCD_LO"  ,"h_QCD_LO FEWZ  ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_QCD_LO  ->Sumw2(); 
   TH1D *h_QCD_NNLO = new TH1D("h_QCD_NNLO","h_QCD_NNLO FEWZ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_QCD_NNLO->Sumw2(); 

   TH1D *h_e_Comb_QCD_EW_Add = new TH1D("h_e_Comb_QCD_EW_Add","h_e_Comb_QCD_EW_Additive  ; #it{M}_{inv} (GeV); d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_Comb_QCD_EW_Add->Sumw2();
   TH1D *h_e_Comb_QCD_EW_Fac = new TH1D("h_e_Comb_QCD_EW_Fac","h_e_Comb_QCD_EW_Factorized; #it{M}_{inv} (GeV); d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_Comb_QCD_EW_Fac->Sumw2();

   TH1D *h_e_EW_Gap   = new TH1D("h_e_EW_Gap"  ,"h_e_EW Gap = NLO - LO (MCSANC)  ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_EW_Gap  ->Sumw2(); 
   TH1D *h_e_EW_Ratio = new TH1D("h_e_EW_Ratio","h_e_EW Ratio = NLO / LO (MCSANC); #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_EW_Ratio->Sumw2(); 

   TH1D *h_e_kfac_add = new TH1D("h_e_kfac_add","K-factor_Additive_ele  ; #it{M}_{inv} (GeV) ; k-factor (/100 GeV)",8000,0,8000);      h_e_kfac_add->Sumw2(); 
   TH1D *h_e_kfac_fac = new TH1D("h_e_kfac_fac","K-factor_Factorized_ele; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_kfac_fac->Sumw2(); 
   TH1D *h_e_kfac_mid = new TH1D("h_e_kfac_mid","K-factor_Mean_ele      ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);  h_e_kfac_mid->Sumw2(); 

   // EWK (MCSANC) array =>  mass(*left edge* value of the bin) xsec(weight)  error
   double p_mi , p_weight , p_error ;
   double m_mi , m_weight , m_error ;
   double nlo_p_mi , nlo_p_weight , nlo_p_error ;
   double nlo_m_mi , nlo_m_weight , nlo_m_error ;
   double ewk_lo_weight;       double ewk_nlo_weight;       
//   double p_mi;         double nlo_p_mi;         
//   double p_weight ;    double nlo_p_weight ;    
//   double p_error ;     double nlo_p_error ;     
//   double m_mi;         double nlo_m_mi;         
//   double m_weight ;    double nlo_m_weight ;    
//   double m_error ;     double nlo_m_error ;     

   // QCD (FEWZ) array => mass(*middle* value of the bin) weight error
   double lo_p_mi , lo_p_weight , lo_p_error ;
   double lo_m_mi , lo_m_weight , lo_m_error ;
   double nnlo_p_mi , nnlo_p_weight , nnlo_p_error ;
   double nnlo_m_mi , nnlo_m_weight , nnlo_m_error ;
   double qcd_lo_weight;    double qcd_nnlo_weight;

//   double lo_p_mi;      double nnlo_p_mi;
//   double lo_p_weight ; double nnlo_p_weight ;
//   double lo_p_error ;  double nnlo_p_error ;
//   double lo_m_mi;      double nnlo_m_mi;
//   double lo_m_weight ; double nnlo_m_weight ;
//   double lo_m_error ;  double nnlo_m_error ;

   double mi; double lo_mi; double nlo_mi; double nnlo_mi ;
   double e_k_mid ;

   int nlines = 0 ;
   int nlinesQCD = 0;
   // TFile* output = new TFile("WToLNu_QCD_EW_StoredHists_NNPDF31_0603.root","RECREATE");

   while(1){
        inp_ewk_lo >> p_mi >> p_weight >> p_error ;
        inm_ewk_lo >> m_mi >> m_weight >> m_error ;
        inp_ewk_nlo >> nlo_p_mi >> nlo_p_weight >> nlo_p_error ;
        inm_ewk_nlo >> nlo_m_mi >> nlo_m_weight >> nlo_m_error ;

	ewk_lo_weight   = p_weight      +  m_weight      ;
	ewk_nlo_weight  = nlo_p_weight  +  nlo_m_weight  ;

	if (!inp_ewk_lo.good()) break;
	h_e_EW_LO  ->SetBinContent(p_mi, ewk_lo_weight);
	h_e_EW_LO  ->SetBinError(p_mi, p_error);
	h_e_EW_NLO ->SetBinContent(p_mi, ewk_nlo_weight);
	h_e_EW_NLO ->SetBinError(p_mi, nlo_p_error);

	nlines++;
   }
   while(1){
	inp_qcd_lo >> lo_p_mi >> lo_p_weight >> lo_p_error ;
	inm_qcd_lo >> lo_m_mi >> lo_m_weight >> lo_m_error ;
	inp_qcd_nnlo >> nnlo_p_mi >> nnlo_p_weight >> nnlo_p_error ;
	inm_qcd_nnlo >> nnlo_m_mi >> nnlo_m_weight >> nnlo_m_error ;

	qcd_lo_weight   = lo_p_weight   +  lo_m_weight   ;
	qcd_nnlo_weight = nnlo_p_weight +  nnlo_m_weight ;

	if (!inp_qcd_lo.good()) break;

        // QCD LO,NNLO FEWZ data have *different mass binning* for three ranges..(due to huge timing in heavy cases)
        // 1GeV bin for mass 50-99 // 5GeV for m 102-997 // 50GeV for m 1025-7975
        // *mi is not the left edge of the bin. it's middle value of the bin.
        // making qcd_mi to the left edge of the bin like MCSANC style. 
	if (lo_p_mi < 100.0) {  
             mi = lo_p_mi - 0.5  ; 
	}else if (lo_p_mi < 1000.0) { 
             mi = lo_p_mi -2.5;
	}else { 
             mi = lo_p_mi - 25.0 ; 
        }
	//cout << "m " << mi << "	lo " << lo_weight << "	nnlo="<< nnlo_weight << endl;
	//if (lo_p_mi != nnlo_p_mi) { cout << "Different Mass " << lo_p_mi << " " << nnlo_p_mi << endl; }
	h_QCD_LO  ->SetBinContent(lo_p_mi, qcd_lo_weight);
	h_QCD_LO  ->SetBinError(lo_p_mi, lo_p_error);
	h_QCD_NNLO->SetBinContent(nnlo_p_mi, qcd_nnlo_weight);
	h_QCD_NNLO->SetBinError(nnlo_p_mi, nnlo_p_error);

	nlinesQCD++;
   }

   h_base_lo->Rebin(50);
   h_base_lo->SetStats(0);
   h_QCD_NNLO->Rebin(100);
   h_QCD_LO  ->Rebin(100);

   h_e_EW_NLO  ->Rebin(100);
   h_e_EW_LO   ->Rebin(100);
   h_e_EW_Gap  ->Rebin(100);
   h_e_EW_Ratio->Rebin(100);
   h_e_Comb_QCD_EW_Add->Rebin(100);
   h_e_Comb_QCD_EW_Fac->Rebin(100);

   h_mu_lo->SetStats(0);
   h_mu_lo  ->Rebin(100);
   h_mu_EW_NLO  ->Rebin(100);
   h_mu_EW_LO   ->Rebin(100);
   h_mu_EW_Gap  ->Rebin(100);
   h_mu_EW_Ratio->Rebin(100);
   h_mu_Comb_QCD_EW_Add->Rebin(100);
   h_mu_Comb_QCD_EW_Fac->Rebin(100);

   h_tau_EW_NLO  ->Rebin(100);
   h_tau_EW_LO   ->Rebin(100);
   h_tau_EW_Gap  ->Rebin(100);
   h_tau_EW_Ratio->Rebin(100);
   h_tau_Comb_QCD_EW_Add->Rebin(100);
   h_tau_Comb_QCD_EW_Fac->Rebin(100);


   TH2D* massframe = new TH2D("massframe", "", 2, 100, 7002, 2, 1E-10, 1E+03);
   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
   massframe->GetXaxis()->SetTitleOffset(1.1);
   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
   massframe->GetYaxis()->SetTitleOffset(1.15);
   massframe->GetXaxis()->SetLabelSize(0.04);
   massframe->GetYaxis()->SetLabelSize(0.04);
   massframe->GetXaxis()->SetTitleSize(0.05);
   massframe->GetYaxis()->SetTitleSize(0.05);
   massframe->SetStats(kFALSE);

   writeExtraText = false; //true;
   extraText = "W#rightarrow#it{e#nu}";
   setTDRStyle();
 
   lumiTextSize =0.5;
   cmsTextSize = 0.5;

        int W = 1200;
        int H = 1200;
        int H_ref = 1200;
        int W_ref = 1200;
        float T = 0.08*H_ref;
        float B = 0.12*H_ref;
        float L = 0.12*W_ref;
        float R = 0.04*W_ref;

//   TCanvas* c1_e = new TCanvas("c1_e","Electron diff. Xsec vs Mlv",800,800);
//   c1_e->SetFillColor(0);
//   c1_e->SetLeftMargin(L/W);
//   c1_e->SetRightMargin(R/W);
//   c1_e->SetTopMargin(T/H);
//   c1_e->SetBottomMargin(B/H);
//   c1_e->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
//   h_base_lo->SetLineWidth(1);      h_e_EW_LO->SetLineColor(kBlack); 
//   h_e_EW_LO->SetLineWidth(1);     h_e_EW_LO->SetLineColor(kGreen); 
//   h_e_EW_NLO->SetLineWidth(1);     h_e_EW_NLO->SetLineColor(kBlue);   
//   h_QCD_LO->SetLineWidth(1);    h_QCD_LO->SetLineColor(kBlack);  
//   h_QCD_NNLO->SetLineWidth(1);     h_QCD_NNLO->SetLineColor(kRed);     
//
//   h_e_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//   h_e_Comb_QCD_EW_Add->SetMarkerSize(1);  
//   h_e_Comb_QCD_EW_Add->SetMarkerStyle(2);  
//   h_e_Comb_QCD_EW_Add->SetMarkerColor(kBlack);  
//
//   h_e_Comb_QCD_EW_Fac->SetMarkerColor(kRed); 
//   h_e_Comb_QCD_EW_Fac->SetMarkerSize(1);   
//   h_e_Comb_QCD_EW_Fac->SetMarkerStyle(3);   
//   h_e_Comb_QCD_EW_Fac->SetMarkerColor(kRed);  
//

   h_e_EW_Gap->Add(h_e_EW_NLO, h_e_EW_LO, 1, -1);
   h_e_EW_Ratio->Divide(h_e_EW_NLO, h_e_EW_LO, 1, 1);

   //h_QCD_LO ->Draw("HIST SAME");
//   h_QCD_NNLO->Draw("HIST SAME");
//   h_e_EW_NLO ->Draw("HIST SAME");
//   h_e_EW_LO  ->Draw("HIST SAME");
//   h_base_lo->Draw("HIST SAME");
//   massframe->Draw("AXISSAME");
//   
//   TLegend* legend = new TLegend(.57,.6,.80,.89, "#it{W}#rightarrow#it{e#nu}");
//   legend->SetTextSize(0.033);
//   legend->SetTextFont(42);
//   legend->AddEntry(h_base_lo     ,"LO Pythia8    "  ,"L"); 
//   legend->AddEntry(h_e_EW_LO  ,"LO EW (MCSANC)"  ,"L"); 
//   legend->AddEntry(h_e_EW_NLO ,"NLO EW (MCSANC)" ,"L"); 
//   //legend->AddEntry(h_QCD_LO ,"LO QCD (FEWZ)"   ,"L"); 
//   legend->AddEntry(h_QCD_NNLO,"NNLO QCD (FEWZ)" ,"L"); 
//   legend->SetBorderSize(0);
//   legend->Draw();
//   
//
//   TCanvas* c1_mu = new TCanvas("c1_mu","Muon diff. Xsec vs Mlv",800,800);
//   c1_mu->SetFillColor(0);
//   c1_mu->SetLeftMargin(L/W);
//   c1_mu->SetRightMargin(R/W);
//   c1_mu->SetTopMargin(T/H);
//   c1_mu->SetBottomMargin(B/H);
//   c1_mu->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
//   //h_mu_lo->SetLineWidth(1);     h_mu_lo->SetLineColor(kBlack); 
//   h_base_lo->SetLineWidth(1);     h_base_lo->SetLineColor(kBlack); 
//   h_mu_EW_LO->SetLineWidth(1);     h_mu_EW_LO->SetLineColor(kGreen); 
//   h_mu_EW_NLO->SetLineWidth(1);     h_mu_EW_NLO->SetLineColor(kBlue);   
//
//   h_mu_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//   h_mu_Comb_QCD_EW_Add->SetMarkerSize(1);  
//   h_mu_Comb_QCD_EW_Add->SetMarkerStyle(2);  
//   h_mu_Comb_QCD_EW_Add->SetMarkerColor(kBlack);  
//
//   h_mu_Comb_QCD_EW_Fac->SetMarkerColor(kRed); 
//   h_mu_Comb_QCD_EW_Fac->SetMarkerSize(1);   
//   h_mu_Comb_QCD_EW_Fac->SetMarkerStyle(3);   
//   h_mu_Comb_QCD_EW_Fac->SetMarkerColor(kRed);  
   h_mu_EW_Gap->Add(h_mu_EW_NLO, h_mu_EW_LO, 1, -1);
   h_mu_EW_Ratio->Divide(h_mu_EW_NLO, h_mu_EW_LO, 1, 1);
//   //h_QCD_LO ->Draw("HIST SAME");
//   h_QCD_NNLO->Draw("HIST SAME");
//   h_mu_EW_NLO ->Draw("HIST SAME");
//   h_mu_EW_LO  ->Draw("HIST SAME");
//   //h_mu_lo->Draw("HIST SAME");
//   h_base_lo->Draw("HIST SAME");
//
//   massframe->Draw("AXISSAME");
//   
//   TLegend* mu_legend = new TLegend(.57,.6,.80,.89, "#it{W}#rightarrow#it{#mu#nu}");
//   mu_legend->SetTextSize(0.033);
//   mu_legend->SetTextFont(42);
//   //mu_legend->AddEntry(h_mu_lo  ,"LO Pythia8"  ,"L"); 
//   mu_legend->AddEntry(h_base_lo  ,"LO Pythia8"  ,"L"); 
//   mu_legend->AddEntry(h_mu_EW_LO  ,"LO EW (MCSANC)"  ,"L"); 
//   mu_legend->AddEntry(h_mu_EW_NLO ,"NLO EW (MCSANC)" ,"L"); 
//   //mu_legend->AddEntry(h_QCD_LO ,"LO QCD (FEWZ)"   ,"L"); 
//   mu_legend->AddEntry(h_QCD_NNLO,"NNLO QCD (FEWZ)" ,"L"); 
//   mu_legend->SetBorderSize(0);
//   mu_legend->Draw();
//
//
//   TCanvas* c1_tau = new TCanvas("c1_tau","Tau diff. Xsec vs Mlv",800,800);
//   c1_tau->SetFillColor(0);
//   c1_tau->SetLeftMargin(L/W);
//   c1_tau->SetRightMargin(R/W);
//   c1_tau->SetTopMargin(T/H);
//   c1_tau->SetBottomMargin(B/H);
//   c1_tau->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
//   h_tau_EW_LO->SetLineWidth(1);     h_tau_EW_LO->SetLineColor(kGreen); 
//   h_tau_EW_NLO->SetLineWidth(1);     h_tau_EW_NLO->SetLineColor(kBlue);   
//
//   h_tau_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//   h_tau_Comb_QCD_EW_Add->SetMarkerSize(1);  
//   h_tau_Comb_QCD_EW_Add->SetMarkerStyle(2);  
//   h_tau_Comb_QCD_EW_Add->SetMarkerColor(kBlack);  
//
//   h_tau_Comb_QCD_EW_Fac->SetMarkerColor(kRed); 
//   h_tau_Comb_QCD_EW_Fac->SetMarkerSize(1);   
//   h_tau_Comb_QCD_EW_Fac->SetMarkerStyle(3);   
//   h_tau_Comb_QCD_EW_Fac->SetMarkerColor(kRed);  

   h_tau_EW_Gap->Add(h_tau_EW_NLO, h_tau_EW_LO, 1, -1);
   h_tau_EW_Ratio->Divide(h_tau_EW_NLO, h_tau_EW_LO, 1, 1);

//   //h_QCD_LO ->Draw("HIST SAME");
//   h_QCD_NNLO->Draw("HIST SAME");
//   h_tau_EW_NLO ->Draw("HIST SAME");
//   h_tau_EW_LO  ->Draw("HIST SAME");
//   h_base_lo->Draw("HIST SAME");
//
//   massframe->Draw("AXISSAME");
//   
//   TLegend* tau_legend = new TLegend(.57,.6,.80,.89, "#it{W}#rightarrow#it{#tau#nu}");
//   tau_legend->SetTextSize(0.033);
//   tau_legend->SetTextFont(42);
//   tau_legend->AddEntry(h_base_lo  ,"LO Pythia8"  ,"L"); 
//   tau_legend->AddEntry(h_tau_EW_LO  ,"LO EW (MCSANC)"  ,"L"); 
//   tau_legend->AddEntry(h_tau_EW_NLO ,"NLO EW (MCSANC)" ,"L"); 
//   //tau_legend->AddEntry(h_QCD_LO ,"LO QCD (FEWZ)"   ,"L"); 
//   tau_legend->AddEntry(h_QCD_NNLO,"NNLO QCD (FEWZ)" ,"L"); 
//   tau_legend->SetBorderSize(0);
//   tau_legend->Draw();
//
//   CMS_lumi(c1_e,14,11);
//   c1_e->Update();
//   c1_e->RedrawAxis();
//   c1_e->GetFrame()->Draw();  c1_e->RedrawAxis();
//   c1_e->GetFrame()->Draw();
//
//   CMS_lumi(c1_mu,14,11);
//   c1_mu->Update();
//   c1_mu->RedrawAxis();
//   c1_mu->GetFrame()->Draw();  c1_e->RedrawAxis();
//   c1_mu->GetFrame()->Draw();
//
//   CMS_lumi(c1_tau,14,11);
//   c1_tau->Update();
//   c1_tau->RedrawAxis();
//   c1_tau->GetFrame()->Draw();  c1_e->RedrawAxis();
//   c1_tau->GetFrame()->Draw();

   //Result
//   TCanvas* c2_e = new TCanvas("c2_e","Electron QCD EW Combined Corrections",800,800);
//   c2_e->SetFillColor(0);
//   c2_e->SetLeftMargin(L/W);
//   c2_e->SetRightMargin(R/W);
//   c2_e->SetTopMargin(T/H);
//   c2_e->SetBottomMargin(B/H);
//   c2_e->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
//   h_base_lo->SetLineWidth(1);     h_base_lo->SetLineColor(kAzure+1); 
   h_e_Comb_QCD_EW_Add->Add(h_e_EW_Gap, h_QCD_NNLO, 1, 1);
   h_e_Comb_QCD_EW_Fac->Multiply(h_e_EW_Ratio, h_QCD_NNLO, 1, 1);

//   h_e_EW_Gap->SetLineColor(kBlue);
//   h_e_EW_Ratio->SetLineColor(kGreen);   
//
//   h_e_Comb_QCD_EW_Fac->SetLineColor(kRed);  
//   h_e_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//
//   //h_e_EW_Gap->Draw("HIST SAME");
//   //h_e_EW_Ratio->Draw("HIST SAME");
//
//   h_e_Comb_QCD_EW_Add->Draw("HIST SAME");
//   h_e_Comb_QCD_EW_Fac->Draw("HIST SAME");
//   h_base_lo->Draw("HIST SAME E");
//   massframe->Draw("AXISSAME");
//
//   TLegend* legend2 = new TLegend(.17,.2,.30,.45, "#it{W}#rightarrow#it{e#nu}");
//   legend2->SetTextSize(0.030);
//   legend2->SetTextFont(42);
//   legend2->AddEntry(h_base_lo, "LO Pythia" ,"LP"); 
//   legend2->AddEntry(h_e_Comb_QCD_EW_Fac, "QCD #otimes EW" ,"L"); 
//   legend2->AddEntry(h_e_Comb_QCD_EW_Add, "QCD #oplus EW" ,"L"); 
//   //legend2->AddEntry(h_e_EW_Gap,  "EW NLO-LO " ,"L"); 
//   //legend2->AddEntry(h_e_EW_Ratio,"EW NLO/LO " ,"L"); 
//   legend2->SetBorderSize(0);
//   legend2->Draw();
//
//   TCanvas* c2_mu = new TCanvas("c2_mu","Muon QCD EW Combined Corrections",800,800);
//   c2_mu->SetFillColor(0);
//   c2_mu->SetLeftMargin(L/W);
//   c2_mu->SetRightMargin(R/W);
//   c2_mu->SetTopMargin(T/H);
//   c2_mu->SetBottomMargin(B/H);
//   c2_mu->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
   h_mu_Comb_QCD_EW_Add->Add(h_mu_EW_Gap, h_QCD_NNLO, 1, 1);
   h_mu_Comb_QCD_EW_Fac->Multiply(h_mu_EW_Ratio, h_QCD_NNLO, 1, 1);

//   //h_mu_lo->SetLineWidth(1);     h_mu_lo->SetLineColor(kAzure+1); 
//   h_base_lo->SetLineWidth(1);     h_base_lo->SetLineColor(kAzure+1); 
//   h_mu_EW_Gap->SetLineColor(kBlue);
//   h_mu_EW_Ratio->SetLineColor(kGreen);   
//
//   h_mu_Comb_QCD_EW_Fac->SetLineColor(kRed);  
//   h_mu_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//
//   //h_mu_EW_Gap->Draw("HIST SAME");
//   //h_mu_EW_Ratio->Draw("HIST SAME");
//
//   h_mu_Comb_QCD_EW_Add->Draw("HIST SAME");
//   h_mu_Comb_QCD_EW_Fac->Draw("HIST SAME");
//   //h_mu_lo->Draw("HIST SAME E");
//   h_base_lo->Draw("HIST SAME E");
//   massframe->Draw("AXISSAME");
//
//   TLegend* mu_legend2 = new TLegend(.17,.2,.30,.45, "#it{W}#rightarrow#it{#mu#nu}");
//   mu_legend2->SetTextSize(0.030);
//   mu_legend2->SetTextFont(42);
//   //mu_legend2->AddEntry(h_mu_lo,"LO Pythia8 " ,"LP"); 
//   mu_legend2->AddEntry(h_base_lo,"LO Pythia8 " ,"LP"); 
//   mu_legend2->AddEntry(h_mu_Comb_QCD_EW_Fac, "QCD #otimes EW" ,"L"); 
//   mu_legend2->AddEntry(h_mu_Comb_QCD_EW_Add, "QCD #oplus EW" ,"L"); 
//   //mu_legend2->AddEntry(h_mu_EW_Gap,  "EW NLO-LO " ,"L"); 
//   //mu_legend2->AddEntry(h_mu_EW_Ratio,"EW NLO/LO " ,"L"); 
//   mu_legend2->SetBorderSize(0);
//   mu_legend2->Draw();
//
//   TCanvas* c2_tau = new TCanvas("c2_tau","Tau QCD EW Combined Corrections",800,800);
//   c2_tau->SetFillColor(0);
//   c2_tau->SetLeftMargin(L/W);
//   c2_tau->SetRightMargin(R/W);
//   c2_tau->SetTopMargin(T/H);
//   c2_tau->SetBottomMargin(B/H);
//   c2_tau->SetLogy();
//   massframe->GetXaxis()->SetTitle("#it{M}_{inv} (GeV)");
//   massframe->GetXaxis()->SetTitleOffset(1.1);
//   massframe->GetYaxis()->SetTitle("d#sigma/dM (pb/100 GeV)");
//   massframe->GetXaxis()->SetMoreLogLabels();
//   massframe->GetXaxis()->SetNoExponent(kTRUE);
//   massframe->Draw();
//
   h_tau_Comb_QCD_EW_Add->Add(h_tau_EW_Gap, h_QCD_NNLO, 1, 1);
   h_tau_Comb_QCD_EW_Fac->Multiply(h_tau_EW_Ratio, h_QCD_NNLO, 1, 1);
//
//   h_base_lo->SetLineWidth(1);     h_base_lo->SetLineColor(kAzure+1); 
//   h_tau_EW_Gap->SetLineColor(kBlue);
//   h_tau_EW_Ratio->SetLineColor(kGreen);   
//
//   h_tau_Comb_QCD_EW_Fac->SetLineColor(kRed);  
//   h_tau_Comb_QCD_EW_Add->SetLineColor(kBlack); 
//
//   //h_tau_EW_Gap->Draw("HIST SAME");
//   //h_tau_EW_Ratio->Draw("HIST SAME");
//
//   h_tau_Comb_QCD_EW_Add->Draw("HIST SAME");
//   h_tau_Comb_QCD_EW_Fac->Draw("HIST SAME");
//   h_base_lo->Draw("HIST SAME E");
//   massframe->Draw("AXISSAME");
//
//   TLegend* tau_legend2 = new TLegend(.17,.2,.30,.45, "#it{W}#rightarrow#it{#tau#nu}");
//   tau_legend2->SetTextSize(0.030);
//   tau_legend2->SetTextFont(42);
//   tau_legend2->AddEntry(h_base_lo,"LO Pythia8 " ,"LP"); 
//   tau_legend2->AddEntry(h_tau_Comb_QCD_EW_Fac, "QCD #otimes EW" ,"L"); 
//   tau_legend2->AddEntry(h_tau_Comb_QCD_EW_Add, "QCD #oplus EW" ,"L"); 
//   //tau_legend2->AddEntry(h_tau_EW_Gap,  "EW NLO-LO " ,"L"); 
//   //tau_legend2->AddEntry(h_tau_EW_Ratio,"EW NLO/LO " ,"L"); 
//   tau_legend2->SetBorderSize(0);
//   tau_legend2->Draw();
//   CMS_lumi(c2_e,14,11);
//   c2_e->Update();
//   c2_e->RedrawAxis();
//   c2_e->GetFrame()->Draw();  c2_e->RedrawAxis();
//   c2_e->GetFrame()->Draw();
//   CMS_lumi(c2_tau,14,11);
//   c2_tau->Update();
//   c2_tau->RedrawAxis();
//   c2_tau->GetFrame()->Draw();  c2_tau->RedrawAxis();
//   c2_tau->GetFrame()->Draw();
//
//   //Write
//   h_QCD_LO  ->Write();
//   h_QCD_NNLO->Write();
//   h_base_lo  ->Write();
//   h_e_EW_NLO  ->Write();
//   h_e_EW_LO   ->Write();
//   h_e_EW_Gap->Write();   
//   h_e_EW_Ratio->Write();  
//   h_e_Comb_QCD_EW_Add->Write();   
//   h_e_Comb_QCD_EW_Fac->Write();  
//
//   c1_e->Write();
//   c2_e->Write();
//
//   //h_mu_lo  ->Write();
//   h_mu_EW_NLO  ->Write();
//   h_mu_EW_LO   ->Write();
//   h_mu_EW_Gap->Write();   
//   h_mu_EW_Ratio->Write();  
//   h_mu_Comb_QCD_EW_Add->Write();   
//   h_mu_Comb_QCD_EW_Fac->Write();  
//
//   c1_mu->Write();
//   c2_mu->Write();
//
//   h_base_lo  ->Write();
//   h_tau_EW_NLO  ->Write();
//   h_tau_EW_LO   ->Write();
//   h_tau_EW_Gap->Write();   
//   h_tau_EW_Ratio->Write();  
//   h_tau_Comb_QCD_EW_Add->Write();   
//   h_tau_Comb_QCD_EW_Fac->Write();  
//
//   c1_tau->Write();
//   c2_tau->Write();
//





  //K-factor
   TH2D* null = new TH2D("null", "K-factor", 2, 101, 7002, 2, 0.0, 1.3 );
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->GetYaxis()->SetTitleOffset(1.15);
   null->GetXaxis()->SetLabelSize(0.04);
   null->GetYaxis()->SetLabelSize(0.04);
   null->GetXaxis()->SetTitleSize(0.05);
   null->GetYaxis()->SetTitleSize(0.05);
   null->SetStats(kFALSE);

   h_e_kfac_add->SetLineColor(kBlack);       h_e_kfac_add->SetLineColor(kBlack);
   h_e_kfac_add->SetMarkerColor(kBlack);     h_e_kfac_add->SetMarkerColor(kBlack);
   h_e_kfac_add->SetMarkerStyle(2);          h_e_kfac_add->SetMarkerStyle(2);    
   h_e_kfac_add->SetMarkerSize(1);           h_e_kfac_add->SetMarkerSize(1);
   h_e_kfac_fac->SetLineColor(kRed);         h_e_kfac_fac->SetLineColor(kRed);
   h_e_kfac_fac->SetMarkerColor(kRed);       h_e_kfac_fac->SetMarkerColor(kRed);  
   h_e_kfac_fac->SetMarkerStyle(3);          h_e_kfac_fac->SetMarkerStyle(3);
   h_e_kfac_fac->SetMarkerSize(1);           h_e_kfac_fac->SetMarkerSize(1);

   h_m_kfac_add->SetLineColor(kBlack);       h_m_kfac_add->SetLineColor(kBlack);
   h_m_kfac_add->SetMarkerColor(kBlack);     h_m_kfac_add->SetMarkerColor(kBlack);
   h_m_kfac_add->SetMarkerStyle(2);          h_m_kfac_add->SetMarkerStyle(2);    
   h_m_kfac_add->SetMarkerSize(1);           h_m_kfac_add->SetMarkerSize(1);
   h_m_kfac_fac->SetLineColor(kRed);         h_m_kfac_fac->SetLineColor(kRed);
   h_m_kfac_fac->SetMarkerColor(kRed);       h_m_kfac_fac->SetMarkerColor(kRed);  
   h_m_kfac_fac->SetMarkerStyle(3);          h_m_kfac_fac->SetMarkerStyle(3);
   h_m_kfac_fac->SetMarkerSize(1);           h_m_kfac_fac->SetMarkerSize(1);

   h_t_kfac_add->SetLineColor(kBlack);       h_t_kfac_add->SetLineColor(kBlack);
   h_t_kfac_add->SetMarkerColor(kBlack);     h_t_kfac_add->SetMarkerColor(kBlack);
   h_t_kfac_add->SetMarkerStyle(2);          h_t_kfac_add->SetMarkerStyle(2);    
   h_t_kfac_add->SetMarkerSize(1);           h_t_kfac_add->SetMarkerSize(1);
   h_t_kfac_fac->SetLineColor(kRed);         h_t_kfac_fac->SetLineColor(kRed);
   h_t_kfac_fac->SetMarkerColor(kRed);       h_t_kfac_fac->SetMarkerColor(kRed);  
   h_t_kfac_fac->SetMarkerStyle(3);          h_t_kfac_fac->SetMarkerStyle(3);
   h_t_kfac_fac->SetMarkerSize(1);           h_t_kfac_fac->SetMarkerSize(1);

   h_e_kfac_add->Rebin(100);   
   h_e_kfac_fac->Rebin(100); 
   h_e_kfac_mid->Rebin(100);   
   h_m_kfac_add->Rebin(100); 
   h_m_kfac_fac->Rebin(100); 
   h_m_kfac_mid->Rebin(100); 
   h_t_kfac_add->Rebin(100); 
   h_t_kfac_fac->Rebin(100); 
   //h_t_kfac_mid->Rebin(100); 

   h_e_kfac_add->Divide(h_e_Comb_QCD_EW_Add,  h_base_lo, 1, 1, "B") ;   
   h_e_kfac_fac->Divide(h_e_Comb_QCD_EW_Fac,  h_base_lo, 1, 1, "B")  ; 

   h_m_kfac_add->Divide(h_mu_Comb_QCD_EW_Add, h_base_lo, 1, 1, "B")  ; 
   h_m_kfac_fac->Divide(h_mu_Comb_QCD_EW_Fac, h_base_lo, 1, 1, "B")  ; 

   h_t_kfac_add->Divide(h_tau_Comb_QCD_EW_Add, h_base_lo, 1, 1, "B")  ; 
   h_t_kfac_fac->Divide(h_tau_Comb_QCD_EW_Fac, h_base_lo, 1, 1, "B")  ; 


   //Fitting example(pol4 fit)
   //TF1 *fitFcn = new TF1("fitFcn", "pol4", 100, 7000); //x range ,4 parameters
   //fitFcn->SetLineColor(kMagenta);
   //h_e_kfac_mid->Add(h_m_kfac_add, h_m_kfac_fac, 1, 1) ;   
   //h_e_kfac_mid->Scale(0.5) ;  
   //h_e_kfac_mid->Draw("same");
   //h_e_kfac_mid->Fit("fitFcn","V+","ep"); 
   //h_e_kfac_mid->Fit(""); 


   TGraph* gr_e_kfac_mid = new TGraph(0);
   TGraphAsymmErrors* err_e_kfac_mid = new TGraphAsymmErrors(0);

   TGraph* gr_m_kfac_mid = new TGraph(0);
   TGraphAsymmErrors* err_m_kfac_mid = new TGraphAsymmErrors(0);

   TGraph* gr_e_kfac_add = new TGraph(0);
   TGraphAsymmErrors* err_e_kfac_add = new TGraphAsymmErrors(0);

   TGraph* gr_m_kfac_add = new TGraph(0);
   TGraphAsymmErrors* err_m_kfac_add = new TGraphAsymmErrors(0);

   cout << "Totla Bin Number " << h_e_kfac_add->GetNbinsX() << endl;

   for(int i=1; i<h_e_kfac_add->GetNbinsX(); i++){

	double maxBin      = ( h_e_kfac_add->GetBinCenter(i) + h_e_kfac_add->GetBinWidth(i) / 2.0 );
        double e_k_mid = getEKfactorMid(maxBin, 0); 
        double m_k_mid = getMKfactorMid(maxBin, 0); 
        double e_k_mid_add = getEKfactorMid(maxBin, 1); 
        double m_k_mid_add = getMKfactorMid(maxBin, 1); 

	double e_k_add     = h_e_kfac_add->GetBinContent(i) ;
	double e_k_add_err = h_e_kfac_add->GetBinError(i) ;
	double e_k_fac     = h_e_kfac_fac->GetBinContent(i) ;
	double e_k_fac_err = h_e_kfac_fac->GetBinError(i) ;

	double m_k_add     = h_m_kfac_add->GetBinContent(i) ;
	double m_k_add_err = h_m_kfac_add->GetBinError(i) ;
	double m_k_fac     = h_m_kfac_fac->GetBinContent(i) ;
	double m_k_fac_err = h_m_kfac_fac->GetBinError(i) ;

	double t_k_add     = h_t_kfac_add->GetBinContent(i) ;
	double t_k_add_err = h_t_kfac_add->GetBinError(i) ;
	double t_k_fac     = h_t_kfac_fac->GetBinContent(i) ;
	double t_k_fac_err = h_t_kfac_fac->GetBinError(i) ;

        gr_e_kfac_mid->SetPoint(i, maxBin, e_k_mid);
        err_e_kfac_mid->SetPoint(i, maxBin, e_k_mid);

        gr_m_kfac_mid->SetPoint(i, maxBin, m_k_mid);
        err_m_kfac_mid->SetPoint(i, maxBin, m_k_mid);

        gr_e_kfac_add->SetPoint(i, maxBin, e_k_mid_add);
        err_e_kfac_add->SetPoint(i, maxBin, e_k_mid_add);

        gr_m_kfac_add->SetPoint(i, maxBin, m_k_mid_add);
        err_m_kfac_add->SetPoint(i, maxBin, m_k_mid_add);

	if(i<2){ 
		err_e_kfac_mid->SetPointError(i, 50., 50., 0.025, 0.025); 
		err_m_kfac_mid->SetPointError(i, 50., 50., 0.025, 0.025); 

		err_e_kfac_add->SetPointError(i, 50., 50., 0.025, 0.025); 
		err_m_kfac_add->SetPointError(i, 50., 50., 0.025, 0.025); 
	}else{
		double e_k_mid_err_max = ( max( abs(e_k_add - e_k_mid) , abs(e_k_fac - e_k_mid) ) > e_k_mid_err_max ) ?  max( abs(e_k_add - e_k_mid) , abs(e_k_fac - e_k_mid) ) : e_k_mid_err_max ;
		double m_k_mid_err_max = ( max( abs(m_k_add - m_k_mid) , abs(m_k_fac - m_k_mid) ) > m_k_mid_err_max ) ?  max( abs(m_k_add - m_k_mid) , abs(m_k_fac - m_k_mid) ) : m_k_mid_err_max ;
		if (e_k_mid_err_max > 0.25) {cout << i << " Big Uncert>0.25 -> e mass " << maxBin << " mid_err " <<  e_k_mid_err_max << endl;  e_k_mid_err_max = 0.24;} 
		if ( i== 36) {m_k_mid_err_max = 0.09;}
		if ( i== 44) {m_k_mid_err_max = 0.18;}
		if (m_k_mid_err_max > 0.19) {cout << i << " Big Uncert>0.19 -> m mass " << maxBin << " mid_err " <<  m_k_mid_err_max << endl;   m_k_mid_err_max = 0.18;}

		err_e_kfac_mid->SetPointError(i, 50., 50., e_k_mid_err_max, e_k_mid_err_max);
		err_m_kfac_mid->SetPointError(i, 50., 50., m_k_mid_err_max, m_k_mid_err_max);

		double e_k_mid_err_max_a = ( abs(e_k_fac - e_k_mid_add) > e_k_mid_err_max_a ) ?  abs(e_k_fac - e_k_mid_add) : e_k_mid_err_max_a ;
		double m_k_mid_err_max_a = ( abs(m_k_fac - m_k_mid_add) > m_k_mid_err_max_a ) ?  abs(m_k_fac - m_k_mid_add) : m_k_mid_err_max_a ;

		if (e_k_mid_err_max_a > 0.25) {cout << i << " Big Uncert>0.25 -> e mass " << maxBin << " mid_err " <<  e_k_mid_err_max_a << endl;  e_k_mid_err_max_a = 0.24;} 
		if ( i== 36) {m_k_mid_err_max_a = 0.09;}
		if ( i== 44) {m_k_mid_err_max_a = 0.19;}
		if (m_k_mid_err_max_a > 0.19) {cout << i << " Big Uncert>0.19 -> m mass " << maxBin << " mid_err " <<  m_k_mid_err_max_a << endl;  m_k_mid_err_max_a = 0.16;}

		err_e_kfac_add->SetPointError(i, 50., 50., e_k_mid_err_max_a, e_k_mid_err_max_a);
		err_m_kfac_add->SetPointError(i, 50., 50., m_k_mid_err_max_a, m_k_mid_err_max_a);
	}
	double e_k_mid_err1 = abs(e_k_add - e_k_mid);
	double e_k_mid_err2 = abs(e_k_fac - e_k_mid);
	double e_k_mid_err_max_a = abs(e_k_fac - e_k_mid);
	double e_k_mid_err_max = max(e_k_mid_err1, e_k_mid_err2);

	double m_k_mid_err1 = abs(m_k_add - m_k_mid);
	double m_k_mid_err2 = abs(m_k_fac - m_k_mid);
	double m_k_mid_err_max_a = abs(m_k_fac - m_k_mid);
	double m_k_mid_err_max = max(m_k_mid_err1, m_k_mid_err2);

	cout << i << " Bin " << maxBin << "********************************************************" << endl;
        cout << "el add=" << e_k_add << "+-" << e_k_add_err << " fac=" << e_k_fac << "+-" << e_k_fac_err << endl;
        cout << "## Fit Mid=" << e_k_mid <<  " (M-add) " << e_k_mid_err1 << " (M-fac) " << e_k_mid_err2 << " Max " << e_k_mid_err_max << " ### e Unc=" << (e_k_mid_err_max / e_k_mid) << endl;
        cout << "## Fit Add=" << e_k_mid_add <<  " (add-fac) " << e_k_mid_err_max_a << " ### e Unc=" << (e_k_mid_err_max_a / e_k_mid_add) << endl;
	cout <<"--------------------------------------------------------------" << endl;
        cout << "mu add=" << m_k_add << "+-" << m_k_add_err << " fac=" << m_k_fac << "+-" << m_k_fac_err << endl;
        cout << "## Fit Mid=" << m_k_mid <<  " (M-add) " << m_k_mid_err1 << " (M-fac) " << m_k_mid_err2 << " Max " << m_k_mid_err_max << " ### m Unc=" << (m_k_mid_err_max / m_k_mid) << endl;
        cout << "## Fit Add=" << m_k_mid_add <<  " (add-fac) " << m_k_mid_err_max_a << " ### m Unc=" << (m_k_mid_err_max_a / m_k_mid_add) << endl;

   }//for binNx

   TFile* outputfile = new TFile("hist_kfactor.root","recreate");

   TCanvas* c3_e_result = new TCanvas("c3_e_result","Electron K-factor (Fit Median) ",800,800);
   c3_e_result->SetFillColor(0);
   c3_e_result->SetLeftMargin(L/W);
   c3_e_result->SetRightMargin(R/W);
   c3_e_result->SetTopMargin(T/H);
   c3_e_result->SetBottomMargin(B/H);
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->Draw();

   gr_e_kfac_mid->SetLineColor(kBlue);
   gr_e_kfac_mid->SetLineWidth(2);

   err_e_kfac_mid->SetFillColorAlpha(kGray, 0.8);
   err_e_kfac_mid->SetFillStyle(3001);
   err_e_kfac_mid->SetLineWidth(0);

   h_e_kfac_add->Draw("SAME E");
   gr_e_kfac_mid->Draw("SAMEC"); 
   err_e_kfac_mid->Draw("SAME E3");
   h_e_kfac_fac->Draw("SAME E");
   null->Draw("AXISSAME");

   TLatex* text = new TLatex();
   text->SetNDC();
   text->SetTextColor(kBlue);
   text->SetTextSize(0.025);
   text->DrawLatex(0.58,0.55,"Fit Parmeters(Median)");   
   text->DrawLatex(0.58,0.50,"Chi2 = 430.183");   
   text->DrawLatex(0.58,0.45,"NDf  = 66");   
   text->DrawLatex(0.58,0.40,"p0   = 1.21213 #pm 6.85406e-04");   
   text->DrawLatex(0.58,0.35,"p1   = -1.9292-e04 #pm 4.44512e-06");   
   text->DrawLatex(0.58,0.30,"p2   = 7.02667e-08 #pm 4.03942e-09");   
   text->DrawLatex(0.58,0.25,"p3   = -1.09576e-11 #pm 1.13511e-12");   
   text->DrawLatex(0.58,0.20,"p4   = 5.17257e-16 #pm 9.60821e-17");   

   TLegend* legend3 = new TLegend(.20,0.2,0.55,0.55, "#it{W}#rightarrow#it{e#nu}");
   legend3->SetTextSize(0.035);
   legend3->SetTextFont(42);
   legend3->AddEntry(h_e_kfac_fac, "QCD #otimes EW" ,"L"); 
   legend3->AddEntry(h_e_kfac_add, "QCD #oplus EW" ,"L"); 
   legend3->AddEntry(gr_e_kfac_mid, "Parameterization" ,"L"); 
   legend3->AddEntry(err_e_kfac_mid, "Uncertainty" ,"F"); 
   legend3->SetBorderSize(0);
   legend3->Draw();

   TCanvas* c3_e_result_add = new TCanvas("c3_e_result_add","Electron K-factor (Fit Additive) ",800,800);
   c3_e_result_add->SetFillColor(0);
   c3_e_result_add->SetLeftMargin(L/W);
   c3_e_result_add->SetRightMargin(R/W);
   c3_e_result_add->SetTopMargin(T/H);
   c3_e_result_add->SetBottomMargin(B/H);
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->Draw();

   gr_e_kfac_add->SetLineColor(kMagenta);
   gr_e_kfac_add->SetLineWidth(2);

   err_e_kfac_add->SetFillColorAlpha(kGray, 0.8);
   err_e_kfac_add->SetFillStyle(3001);
   err_e_kfac_add->SetLineWidth(0);

   err_e_kfac_add->Draw("SAME E3");
   h_e_kfac_add->Draw("SAME E");
   h_e_kfac_fac->Draw("SAME E");
   gr_e_kfac_add->Draw("SAMEC");
   null->Draw("AXISSAME");

   text->SetTextColor(kMagenta);
   text->SetTextSize(0.025);
   text->DrawLatex(0.58,0.55,"Fit Parameters(Additive)");   
   text->DrawLatex(0.58,0.50,"Chi2 = 451.529");   
   text->DrawLatex(0.58,0.45,"NDf  = 66");   
   text->DrawLatex(0.58,0.40,"p0   = 1.24039 #pm 8.08117e-04");   
   text->DrawLatex(0.58,0.35,"p1   = -2.35493e-04 #pm 4.99893e-06");   
   text->DrawLatex(0.58,0.30,"p2   = 1.03222e-07 #pm 4.10452e-09");   
   text->DrawLatex(0.58,0.25,"p3   = -1.93752e-11 #pm 1.03441e-12");   
   text->DrawLatex(0.58,0.20,"p4   = 1.18182e-15 #pm 7.84593e-17");   

   TLegend* legend31 = new TLegend(.20,0.2,0.55,0.55, "#it{W}#rightarrow#it{e#nu}");
   legend31->SetTextSize(0.035);
   legend31->SetTextFont(42);
   legend31->AddEntry(h_e_kfac_fac, "QCD #otimes EW" ,"L"); 
   legend31->AddEntry(h_e_kfac_add, "QCD #oplus EW" ,"L"); 
   legend31->AddEntry(gr_e_kfac_add, "Parameterization" ,"L"); 
   legend31->AddEntry(err_e_kfac_add, "Uncertainty" ,"F"); 
   legend31->SetBorderSize(0);
   legend31->Draw();

// Fitting example(pol4 fit)
//   TF1 *fitFcn = new TF1("fitFcn", "pol4", 100, 700); //x range ,4 parameters
//   fitFcn->SetLineColor(kMagenta);
//   h_m_kfac_mid->Add(h_m_kfac_add, h_m_kfac_fac, 1, 1) ;   
//   h_m_kfac_mid->Scale(0.5) ;  
//   h_m_kfac_mid->Fit("fitFcn","V+","ep"); 
//   h_e_kfac_mid->Fit(""); 

   TCanvas* c3_m_result = new TCanvas("c3_m_result","Muon K-factor (Fit Median) ",800,800);
   c3_m_result->SetFillColor(0);
   c3_m_result->SetLeftMargin(L/W);
   c3_m_result->SetRightMargin(R/W);
   c3_m_result->SetTopMargin(T/H);
   c3_m_result->SetBottomMargin(B/H);
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->Draw();

   gr_m_kfac_mid->SetLineColor(kBlue);
   gr_m_kfac_mid->SetLineWidth(2);

   err_m_kfac_mid->SetFillColorAlpha(kGray, 0.8);
   err_m_kfac_mid->SetFillStyle(3001);
   err_m_kfac_mid->SetLineWidth(0);

   err_m_kfac_mid->Draw("SAME E3");
   h_m_kfac_add->Draw("SAME E");
   h_m_kfac_fac->Draw("SAME E");
   gr_m_kfac_mid->Draw("SAMEC");
   null->Draw("AXISSAME");


   h_e_kfac_add->Write();
   gr_e_kfac_mid->Write(); 
   h_m_kfac_add->Write(); 
   gr_m_kfac_mid->Write();

   outputfile->Write();
   outputfile->Close();


   text->SetTextColor(kBlue);
   text->SetTextSize(0.025);
   text->DrawLatex(0.58,0.55,"Fit Parmeters(Median)");   
   text->DrawLatex(0.58,0.50,"Chi2 = 322.055");   
   text->DrawLatex(0.58,0.45,"NDf  = 66");   
   text->DrawLatex(0.58,0.40,"p0   = 1.21384 #pm 7.20025e-04");   
   text->DrawLatex(0.58,0.35,"p1   = -2.0136e-04 #pm 4.75479e-06");   
   text->DrawLatex(0.58,0.30,"p2   = 8.32023e-08 #pm 4.43062e-09");   
   text->DrawLatex(0.58,0.25,"p3   = -1.55027e-11 #pm 1.25214e-12");   
   text->DrawLatex(0.58,0.20,"p4   = 9.41259e-16 #pm 1.06378e-16");   

   TLegend* legendm3 = new TLegend(.20,0.2,0.55,0.55, "#it{W}#rightarrow#it{#mu#nu}");
   legendm3->SetTextSize(0.035);
   legendm3->SetTextFont(42);
   legendm3->AddEntry(h_m_kfac_fac, "QCD #otimes EW" ,"L"); 
   legendm3->AddEntry(h_m_kfac_add, "QCD #oplus EW" ,"L"); 
   legendm3->AddEntry(gr_m_kfac_mid, "Parameterization" ,"L"); 
   legendm3->AddEntry(err_m_kfac_mid, "Uncertainty" ,"F"); 
   legendm3->SetBorderSize(0);
   legendm3->Draw();

   TCanvas* c3_m_result_add = new TCanvas("c3_m_result_add","Muon K-factor (Fit Additive) ",800,800);
   c3_m_result_add->SetFillColor(0);
   c3_m_result_add->SetLeftMargin(L/W);
   c3_m_result_add->SetRightMargin(R/W);
   c3_m_result_add->SetTopMargin(T/H);
   c3_m_result_add->SetBottomMargin(B/H);
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->Draw();

   gr_m_kfac_add->SetLineColor(kMagenta);
   gr_m_kfac_add->SetLineWidth(2);

   err_m_kfac_add->SetFillColorAlpha(kGray, 0.8);
   err_m_kfac_add->SetFillStyle(3001);
   err_m_kfac_add->SetLineWidth(0);

   err_m_kfac_add->Draw("SAME E3");
   h_m_kfac_add->Draw("SAME E");
   h_m_kfac_fac->Draw("SAME E");
   gr_m_kfac_add->Draw("SAMEC");
   null->Draw("AXISSAME");

   text->SetTextColor(kMagenta);
   text->SetTextSize(0.025);
   text->DrawLatex(0.58,0.55,"Fit Parameters(Additive)");   
   text->DrawLatex(0.58,0.50,"Chi2 = 204.906");   
   text->DrawLatex(0.58,0.45,"NDf  = 66");   
   text->DrawLatex(0.58,0.40,"p0   = 1.24197 #pm 9.27194e-04");   
   text->DrawLatex(0.58,0.35,"p1   = -2.43528e-04 #pm 6.08478e-06");   
   text->DrawLatex(0.58,0.30,"p2   = 1.15169e-07 #pm 5.65003e-09");   
   text->DrawLatex(0.58,0.25,"p3   = -2.35027e-11 #pm 1.6021e-12");   
   text->DrawLatex(0.58,0.20,"p4   = 1.56328e-15 #pm 1.37149e-16");   


   TLegend* legendm31 = new TLegend(.20,0.2,0.55,0.55, "#it{W}#rightarrow#it{#mu#nu}");
   legendm31->SetTextSize(0.035);
   legendm31->SetTextFont(42);
   legendm31->AddEntry(h_m_kfac_fac, "QCD #otimes EW" ,"L"); 
   legendm31->AddEntry(h_m_kfac_add, "QCD #oplus EW" ,"L"); 
   legendm31->AddEntry(gr_m_kfac_add, "Parameterization" ,"L"); 
   legendm31->AddEntry(err_m_kfac_add, "Uncertainty" ,"F"); 
   legendm31->SetBorderSize(0);
   legendm31->Draw();

//   TCanvas* c3_tau_result = new TCanvas("c3_tau_result","Tau K-factor",800,800);
//   c3_tau_result->SetFillColor(0);
//   c3_tau_result->SetLeftMargin(L/W);
//   c3_tau_result->SetRightMargin(R/W);
//   c3_tau_result->SetTopMargin(T/H);
//   c3_tau_result->SetBottomMargin(B/H);
//   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
//   null->GetXaxis()->SetTitleOffset(1.1);
//   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
//   null->Draw();
//
//   h_t_kfac_add->Draw("SAME E");
//   h_t_kfac_fac->Draw("SAME E");
//   null->Draw("AXISSAME");
//
//   TLegend* tau_legend3 = new TLegend(.2,0.2,0.4,0.4, "#it{W}#rightarrow#it{#tau#nu}");
//   tau_legend3->SetTextSize(0.040);
//   tau_legend3->SetTextFont(42);
//   tau_legend3->AddEntry(h_m_kfac_fac, "QCD #otimes EW" ,"L"); 
//   tau_legend3->AddEntry(h_m_kfac_add, "QCD #oplus EW" ,"L"); 
//   tau_legend3->SetBorderSize(0);
//   tau_legend3->Draw();
//
   CMS_lumi(c3_e_result_add,14,11);
   c3_e_result_add->Update();
   c3_e_result_add->RedrawAxis();
   c3_e_result_add->GetFrame()->Draw();  c3_e_result_add->RedrawAxis();
   c3_e_result_add->GetFrame()->Draw();

   CMS_lumi(c3_e_result,14,11);
   c3_e_result->Update();
   c3_e_result->RedrawAxis();
   c3_e_result->GetFrame()->Draw();  c3_e_result->RedrawAxis();
   c3_e_result->GetFrame()->Draw();

   CMS_lumi(c3_m_result_add,14,11);
   c3_m_result_add->Update();
   c3_m_result_add->RedrawAxis();
   c3_m_result_add->GetFrame()->Draw();  c3_m_result_add->RedrawAxis();
   c3_m_result_add->GetFrame()->Draw();

   CMS_lumi(c3_m_result,14,11);
   c3_m_result->Update();
   c3_m_result->RedrawAxis();
   c3_m_result->GetFrame()->Draw();  c3_m_result->RedrawAxis();
   c3_m_result->GetFrame()->Draw();

//   CMS_lumi(c3_tau_result,14,11);
//   c3_tau_result->Update();
//   c3_tau_result->RedrawAxis();
//   c3_tau_result->GetFrame()->Draw();  c3_tau_result->RedrawAxis();
//   c3_tau_result->GetFrame()->Draw();

//   h_e_kfac_add  ->Write();
//   h_e_kfac_fac  ->Write();
//   c3_e_result->Write();
// 
//   h_m_kfac_add  ->Write();
//   h_m_kfac_fac  ->Write();
//   c3_mu_result->Write();
//
//   h_t_kfac_add  ->Write();
//   h_t_kfac_fac  ->Write();
//   c3_tau_result->Write();
//   output->Print();
//   output->Close();
}
