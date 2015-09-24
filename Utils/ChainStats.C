#include "TMath.h"

void ChainStats(TString ifile = "MCout", int nchains)
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  TCanvas * c = new TCanvas("c","ChainStats",1000,1000);
  //  c->SetFillColor(18);
  //  c->Divide(2,3);
  gStyle->SetOptStat("");

  TFile file(ifile+".root","read");
  TNtupleD *nt0 = file.Get("AnalysisTree");
  TH1D *tmp = new TH1D("tmp","tmp",50, 0., 1000.);

  nt0->Project("tmp","fNParameters");
  long int npars = (long int) tmp->GetMean();
  double ss = TMath::Sqrt(tmp->GetMean());
  int nx = TMath::FloorNint(ss);
  int ny = TMath::CeilNint((double) npars/nx);

  c->Divide(nx,ny);
  
  TNtupleD * nt;
  TGraphErrors tge[npars];
  TGraph tg[npars*nchains];
  Double_t x;

  TCanvas * c2 = new TCanvas("c2","TracePlot",1000,1000);
  c2->Divide(TMath::FloorNint(TMath::Sqrt(nchains)),TMath::CeilNint((double) nchains/TMath::FloorNint(TMath::Sqrt(nchains))));


  for (Int_t i = 0; i < npars; i++)
    {
      c2->SetTitle(Form("TracePlot%d",i));
      for (Int_t j = 0; j < nchains; j++)
	{
	  TString ntp = Form("MarkovChainTree_%d",j);
	  nt = (TNtupleD * ) file.Get(ntp);
	  TString npar = Form("Parameter%d",i);
	  TH1D * histo = new TH1D("histo","",100,nt->GetMinimum(npar),nt->GetMaximum(npar));
	  nt->Project("histo",npar);
	  tge[i].SetPoint(j,j+1,histo->GetMean());
	  tge[i].SetPointError(j,0.,histo->GetRMS());
	  histo->Delete();
	  nt->SetBranchAddress(npar,&x);
	  for(Int_t k = 0; k < nt->GetEntries(); k++){
	    nt->GetEntry(k);
	    tg[i*nchains+j].SetPoint(k,k,x);
	  }
	  tg[i*nchains+j].SetTitle(Form("Chain%d",j));
	  c2->cd(j+1);
	  tg[i*nchains+j].Draw("al");
	}
      c2->cd();
      switch(i){
      case 0:
	c2->Print(ifile+"_traceplots.pdf(","pdf");
	break;
      case (npars-1):
	c2->Print(ifile+"_traceplots.pdf)","pdf");
	break;
      default:
	c2->Print(ifile+"_traceplots.pdf","pdf");
      }
      tge[i].SetTitle(Form("Parameter%d",i));
      c->cd(i+1);
      tge[i].Draw("ap");
    }

  c->cd();
  c->SaveAs(ifile+"_chains.pdf");

  file.Close();
}
