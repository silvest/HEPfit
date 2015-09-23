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
  
  for (Int_t i = 0; i < npars; i++)
    {
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
	}
      tge[i].SetTitle(Form("Parameter%d",i));
      c->cd(i+1);
      tge[i].Draw("ap");
    }

  c->cd();
  c->SaveAs(ifile+"_chains.pdf");

  file.Close();
}
