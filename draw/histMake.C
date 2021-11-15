#include<Riostream.h>
#include "tdrstyle.C"
#include "CMS_lumi.C"

Double_t getEKfactorMid(double x, bool IsAdditive){
	if(x > 7000.0) { x = 7000 ; }
	//Additive Fit
	if(IsAdditive == 1) return (1.24039  -(0.000235493*x) + (1.03222e-07*pow(x,2)) -(1.93752e-11 *pow(x,3)) + (1.18182e-15*pow(x,4)) );
	//Middle of Add Fac Fit
	if(IsAdditive == 0) return ( 1.21213 -(0.00019292*x) + (7.02667e-08*pow(x,2)) -(1.09576e-11*pow(x,3)) + (5.17257e-16*pow(x,4)) );
}

void histMake(){

   TFile* BaseLO   = new TFile("Wmass.root" , "READ");       //from Sebastian
   TH1D* h_base_lo = (TH1D*)BaseLO-> Get("Wgenmass");  

   TString datap_ewk_lo   = "./data/LO_EW_Wp_e_m34_hist_new.dat" ;  //MCSANC LO W+ 
   TString datam_ewk_lo   = "./data/LO_EW_Wm_e_m34_hist_new.dat" ;  //MCSANC LO W-
   TString datap_ewk_nlo  = "./data/NLO_EW_Wp_e_m34_hist_new.dat";  //MCSANC NLO W+ 
   TString datam_ewk_nlo  = "./data/NLO_EW_Wm_e_m34_hist_new.dat";  //MCSANC NLO W- 
   TString datap_qcd_lo   = "./data/LO_QCD_Wp_hist.dat"          ;  //FEWZ LO W+         
   TString datam_qcd_lo   = "./data/LO_QCD_Wm_hist.dat"          ;  //FEWZ LO W-         
   TString datap_qcd_nnlo = "./data/NNLO_QCD_Wp_hist.dat"        ;  //FEWZ NNLO W+       
   TString datam_qcd_nnlo = "./data/NNLO_QCD_Wm_hist.dat"        ;  //FEWZ NNLO W-       

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

   // Define Histograms (all mass bin is 1GeV)
   TH1D *h_e_EW_LO = new TH1D("h_e_EW_LO" ,"h_e_EW_LO MCSANC ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000); h_e_EW_LO ->Sumw2(); 
   TH1D *h_e_EW_NLO= new TH1D("h_e_EW_NLO","h_e_EW_NLO MCSANC; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000); h_e_EW_NLO->Sumw2();
   TH1D *h_QCD_LO  = new TH1D("h_QCD_LO"  ,"h_QCD_LO FEWZ  ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);   h_QCD_LO  ->Sumw2(); 
   TH1D *h_QCD_NNLO= new TH1D("h_QCD_NNLO","h_QCD_NNLO FEWZ; #it{M}_{inv} (GeV) ; d#sigma/dM (pb/100 GeV)",8000,0,8000);   h_QCD_NNLO->Sumw2(); 

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

   // QCD (FEWZ) array => mass(*middle* value of the bin) weight error
   double lo_p_mi , lo_p_weight , lo_p_error ;
   double lo_m_mi , lo_m_weight , lo_m_error ;
   double nnlo_p_mi , nnlo_p_weight , nnlo_p_error ;
   double nnlo_m_mi , nnlo_m_weight , nnlo_m_error ;
   double qcd_lo_weight;    double qcd_nnlo_weight;

   double mi; double lo_mi; double nlo_mi; double nnlo_mi ;
   double e_k_mid ;

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
             mi = lo_p_mi - 2.5;
	}else { 
             mi = lo_p_mi - 25.0 ; 
        }
	//cout << "m " << mi << "	lo " << lo_weight << "	nnlo="<< nnlo_weight << endl;
	//if (lo_p_mi != nnlo_p_mi) { cout << "Different Mass " << lo_p_mi << " " << nnlo_p_mi << endl; }
	h_QCD_LO  ->SetBinContent(lo_p_mi, qcd_lo_weight);
	h_QCD_LO  ->SetBinError(lo_p_mi, lo_p_error);
	h_QCD_NNLO->SetBinContent(nnlo_p_mi, qcd_nnlo_weight);
	h_QCD_NNLO->SetBinError(nnlo_p_mi, nnlo_p_error);
   }

   // Rebinning
   // To make all histograms has the same mass binning of 100 GeV
   h_base_lo   ->Rebin(50);
   h_e_EW_NLO  ->Rebin(100);
   h_e_EW_LO   ->Rebin(100);
   h_e_EW_Gap  ->Rebin(100);
   h_e_EW_Ratio->Rebin(100);
   h_QCD_NNLO->Rebin(100);
   h_QCD_LO  ->Rebin(100);
   h_e_Comb_QCD_EW_Add->Rebin(100);
   h_e_Comb_QCD_EW_Fac->Rebin(100);
   h_e_kfac_add->Rebin(100);   
   h_e_kfac_fac->Rebin(100); 
   h_e_kfac_mid->Rebin(100);   

   h_e_EW_Gap->Add(h_e_EW_NLO, h_e_EW_LO, 1, -1);
   h_e_EW_Ratio->Divide(h_e_EW_NLO, h_e_EW_LO, 1, 1);

   // high order weight of combined QCD+EWK by Add method or Fac method
   h_e_Comb_QCD_EW_Add->Add(h_e_EW_Gap, h_QCD_NNLO, 1, 1);
   h_e_Comb_QCD_EW_Fac->Multiply(h_e_EW_Ratio, h_QCD_NNLO, 1, 1);

   // k-factor with combined QCD and EWK by Add method or Fac method
   h_e_kfac_add->Divide(h_e_Comb_QCD_EW_Add,  h_base_lo, 1, 1, "B") ;   
   h_e_kfac_fac->Divide(h_e_Comb_QCD_EW_Fac,  h_base_lo, 1, 1, "B")  ; 


   TGraph* gr_kfac_fit_add = new TGraph(0);   TGraphAsymmErrors* err_kfac_fit_add = new TGraphAsymmErrors(0);

   for(int i=1; i < h_e_kfac_add->GetNbinsX(); i++){

	double maxBin      = ( h_e_kfac_add->GetBinCenter(i) + h_e_kfac_add->GetBinWidth(i) / 2.0 );
        double e_k_fit_add = getEKfactorMid(maxBin, 1); //parameter fit for additive k-factor 
        double e_k_mid     = getEKfactorMid(maxBin, 0); //parameter fit for mid of add and fac k-factor (not used)

	double e_k_add     = h_e_kfac_add->GetBinContent(i) ;
	double e_k_add_err = h_e_kfac_add->GetBinError(i) ;
	double e_k_fac     = h_e_kfac_fac->GetBinContent(i) ;
	double e_k_fac_err = h_e_kfac_fac->GetBinError(i) ;

        gr_kfac_fit_add->SetPoint(i, maxBin, e_k_fit_add);
        err_kfac_fit_add->SetPoint(i, maxBin, e_k_fit_add);
        double e_k_fit_err_max_a = 0.0;

	if(i<2){ 
		err_kfac_fit_add->SetPointError(i, 50., 50., 0.025, 0.025); 
	}else{
		double e_k_fit_err_max_a =  abs(e_k_fac - e_k_fit_add) ;
		//double e_k_fit_err_max_a = ( abs(e_k_fac - e_k_fit_add) > e_k_fit_err_max_a ) ?  abs(e_k_fac - e_k_fit_add) : e_k_fit_err_max_a ;
                err_kfac_fit_add->SetPointError(i, 50., 50., e_k_fit_err_max_a, e_k_fit_err_max_a);

	        cout << "***** " << i << " bin " << maxBin << "********" << endl;
                //cout << "k-factor_add = " << e_k_add << " +- " << e_k_add_err << " ,   k-factor fac = " << e_k_fac << " +- " << e_k_fac_err << ", diff = " << abs(e_k_fac - e_k_add) <<endl;
                //cout << "paramfit_add = " << e_k_fit_add <<  " ,  paramfit_mid = " << e_k_mid << endl;
                //cout << "diff (addfit - fac) = " << abs(e_k_fac - e_k_fit_add)  << ", (midfit -fac) = " << abs(e_k_fac - e_k_mid) << endl;
                cout << "paramfit_add_err = " << e_k_fit_err_max_a << ",  sys uncert = " << (e_k_fit_err_max_a / e_k_fit_add) * 100 << " %" << endl;
	}
   }//end forloop binNx

   TFile* outputfile = new TFile("hist_kfactor.root","recreate");

   // Drawing K-factor plot
   // setTDRStyle setting
   writeExtraText = false; //true;
   extraText = "W#rightarrow#it{e#nu}";
   setTDRStyle();
   lumiTextSize =0.5;
   cmsTextSize = 0.5;
   int H = 1200;   int H_ref = 1200;
   int W = 1200;   int W_ref = 1200;   
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;

   TCanvas* c3_e_result_add = new TCanvas("c3_e_result_add","Electron K-factor (Fit Additive) ",800,800);
   c3_e_result_add->SetFillColor(0);
   c3_e_result_add->SetLeftMargin(L/W);
   c3_e_result_add->SetRightMargin(R/W);
   c3_e_result_add->SetTopMargin(T/H);
   c3_e_result_add->SetBottomMargin(B/H);
   TH2D* null = new TH2D("null", "K-factor", 2, 101, 7002, 2, 0.0, 1.3 );
   null->GetXaxis()->SetTitle("#it{M}_{W} (GeV)");
   null->GetYaxis()->SetTitle("#sigma((N)NLO)/#sigma(LO)");
   null->GetXaxis()->SetTitleOffset(1.1);
   null->Draw();

   h_e_kfac_add->SetLineColor(kBlack);  h_e_kfac_add->SetMarkerColor(kBlack);     
   h_e_kfac_add->SetMarkerStyle(2);     h_e_kfac_add->SetMarkerSize(1);           

   h_e_kfac_fac->SetLineColor(kRed);    h_e_kfac_fac->SetMarkerColor(kRed);      
   h_e_kfac_fac->SetMarkerStyle(3);     h_e_kfac_fac->SetMarkerSize(1);          

   gr_kfac_fit_add->SetLineColor(kMagenta); 
   gr_kfac_fit_add->SetLineWidth(2);

   err_kfac_fit_add->SetFillColorAlpha(kGray, 0.8);
   err_kfac_fit_add->SetFillStyle(3001);
   err_kfac_fit_add->SetLineWidth(0);

   err_kfac_fit_add->Draw("SAME E3");
   h_e_kfac_add->Draw("SAME E");
   h_e_kfac_fac->Draw("SAME E");
   gr_kfac_fit_add->Draw("SAMEC");
   null->Draw("AXISSAME");
 
   TLatex* text = new TLatex();
   text->SetNDC();
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
   legend31->AddEntry(gr_kfac_fit_add, "Parameterization" ,"L"); 
   legend31->AddEntry(err_kfac_fit_add, "Uncertainty" ,"F"); 
   legend31->SetBorderSize(0);
   legend31->Draw();

   CMS_lumi(c3_e_result_add,14,11);
   c3_e_result_add->Update();
   //c3_e_result_add->RedrawAxis();
   //c3_e_result_add->GetFrame()->Draw(); 


   //Fitting example(pol4 fit) (add+fac/2) case not used
   //TF1 *fitFcn = new TF1("fitFcn", "pol4", 100, 7000); //x range ,4 parameters
   //fitFcn->SetLineColor(kMagenta);
   //h_e_kfac_mid->Add(h_e_kfac_add, h_e_kfac_fac, 1, 1) ;   
   //h_e_kfac_mid->Scale(0.5) ;  
   //h_e_kfac_mid->Draw("same");
   //h_e_kfac_mid->Fit("fitFcn","V+","ep"); 
   //h_e_kfac_mid->Fit(""); 

   //h_e_kfac_add->Write();
   //h_e_kfac_fac->Write();
   //outputfile->Write();
   //outputfile->Close();


}
