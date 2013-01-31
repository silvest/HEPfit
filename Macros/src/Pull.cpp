/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TGaxis.h>
#include <TColor.h>
#include <TLatex.h>
#include "Pull.h"


Pull::Pull(TH1D& hist, const int nbinX, const int nbinY, 
           const double xLow, const double xUp, 
           const double yLow, const double yUp) 
{
    nx = nbinX;
    ny = nbinY;
    x_low = xLow;
    x_up = xUp;
    y_low = yLow;
    y_up = yUp;

    copyHist = (TH1D*) hist.Clone("copy"); 
    double sum = copyHist->Integral();
    copyHist->Scale(1.0/sum);
    
    CompatPlot = NULL;
    lx = new TLine();
    tText = new TLatex();
}


void Pull::Draw(const TString xlab, const TString ylab, 
                const double xval, const double xerr, const int maxDigits, 
                const double YTitleOffset) 
{

    makeCompatPlot();
    
    CompatPlot->GetXaxis()->SetTitleSize(0.06);
    CompatPlot->GetYaxis()->SetTitleSize(0.06);
    CompatPlot->GetXaxis()->SetTitleOffset(1.1);
    CompatPlot->GetYaxis()->SetTitleOffset(YTitleOffset);
    //CompatPlot->GetXaxis()->SetNdivisions(508);
    CompatPlot->GetXaxis()->SetNdivisions(506);
    //CompatPlot->GetYaxis()->SetNdivisions(510);
    CompatPlot->GetYaxis()->SetNdivisions(508);
    CompatPlot->GetXaxis()->SetLabelSize(0.043);
    CompatPlot->GetXaxis()->SetLabelOffset(0.013);
    CompatPlot->GetYaxis()->SetLabelSize(0.043);
    CompatPlot->GetYaxis()->SetLabelOffset(0.013);
    CompatPlot->GetZaxis()->SetLabelSize(0.043);
    CompatPlot->GetZaxis()->SetLabelOffset(0.010);
    CompatPlot->SetLabelFont(42,"X");
    CompatPlot->SetLabelFont(42,"Y");
    CompatPlot->SetLabelFont(42,"Z");
    CompatPlot->SetTitleFont(42,"X");
    CompatPlot->SetTitleFont(42,"Y");
    CompatPlot->SetTitleFont(42,"Z");
    //CompatPlot->SetLabelFont(62,"X");
    //CompatPlot->SetLabelFont(62,"Y");
    //CompatPlot->SetLabelFont(62,"Z");
    //CompatPlot->SetTitleFont(62,"X");
    //CompatPlot->SetTitleFont(62,"Y");
    //CompatPlot->SetTitleFont(62,"Z");
    ((TGaxis*) CompatPlot->GetXaxis())->SetMaxDigits(maxDigits);
    ((TGaxis*) CompatPlot->GetYaxis())->SetMaxDigits(maxDigits);
    
    // Titles of the axes 
    if (xlab.CompareTo("")==0) {
        TString Xtitle = ConvertTitle(copyHist->GetXaxis()->GetTitle());
        CompatPlot->GetXaxis()->SetTitle(Xtitle);
    } else
        CompatPlot->GetXaxis()->SetTitle(xlab);
    if (ylab.CompareTo("")==0) {
        TString Xtitle = CompatPlot->GetXaxis()->GetTitle();
        CompatPlot->GetYaxis()->SetTitle("#sigma(" + Xtitle + ")");
    } else 
        CompatPlot->GetYaxis()->SetTitle(ylab);
    
    // color palette
    TColor col85(65,0.00,0.75,1.00);
    TColor col84(64,0.12,0.56,1.00);
    int colors[6]={0,7,65,64,4,2};
    CompatPlot->SetMaximum(6.0);
    CompatPlot->SetContour(6);
    gStyle->SetPalette(6,colors);
    
    CompatPlot->SetTitle("");
    gPad->SetGrid();
    gROOT->ForceStyle();
    
    CompatPlot->Draw("COLZ");    
    
    // add "sigma" to the label of the z axis
    gPad->SetRightMargin(0.20);
    tText->SetTextSize(0.06);
    double xadd = (x_up - x_low)*0.15;
    double yadd = (y_up - y_low)*0.01;
    tText->DrawLatex(x_up + xadd, y_up - yadd, "#sigma");
    
    // the measured point
    if (xval != -999.0) {
        double lw = (y_up - y_low)/40.0;   
        lx->SetLineWidth(4);
        //std::cout << xval - 0.1 * (y_up - y_low) << std::endl;
        lx->DrawLine(xval - 0.025 * (x_up - x_low), xerr, 
                     xval + 0.025 * (x_up - x_low), xerr);
        lx->DrawLine(xval, xerr + 0.025 * (y_up - y_low), 
                     xval, xerr - 0.025 * (y_up - y_low));
    }

    gPad->Modified();
    gPad->Update();    
}


//////////////////////////////////////////////////////////////////////// 

double Pull::fhisto(const double x) const 
{
    double val;
    int ibin = copyHist->FindBin(x);
    if (ibin < 0 || ibin >= copyHist->GetNbinsX())
        val = 0.0; 
    else
        val = copyHist->GetBinContent(ibin);
    return val;
}


double Pull::fconv(const double x) const 
{
    double f1 = fhisto(x);
    double f2 = TMath::Gaus(x - deltaTmp, meanTmp, sigmaTmp, true);
    if (f1 <= 0.0) f1 = 0.0;
    if (f2 <= 0.0) f2 = 0.0;
    return f1*f2;
}


double Pull::integral(const char* funcname, const double xmin, 
                      const double xmax) 
{
    int n = 100;
    double val = 0.0;
    double x0, x1, x2, w = (xmax - xmin)/(double)n;
    for (int i=0; i<n; i+=2) {
        x0 = xmin + (double)i*w;
        x1 = x0 + w;
        x2 = x1 + w;
        if (strcmp(funcname,"fconv")==0)
            val += w/3.0*(fconv(x0) + 4.0*fconv(x1) + fconv(x2));
        else if (strcmp(funcname,"fdelta")==0)
            val += w/3.0*(fdelta(x0) + 4.0*fdelta(x1) + fdelta(x2));
        else 
            std::cout << "Error in Pull::integral()" << std::endl;
    }
    return val;
}


double Pull::fdelta(const double x) 
{
    deltaTmp = x;
    return integral("fconv", copyHist->GetXaxis()->GetXmin(), 
                    copyHist->GetXaxis()->GetXmax());
}


double Pull::fzero(const double x, const double area) const 
{
    return ( 2.0*area - TMath::Erfc(x/sqrt(2.0)) );
}


double Pull::finverfc(const double x) const 
{
    double w = 3.0, eps = 0.001;
    double diff, val, s = 0.0;

    if (x <= 1.0e-100) return 20.0;
    if (x >= 0.5) return 0.0;
    do {
        if (fzero(s, x) == 0.0) val = s;
        if (fzero(s + w, x) == 0.0) val = s + w;
        diff = fzero(s, x) * fzero(s + w, x);
        if (diff < 0.0) {
            w /= 2.0;
            val = s + w;
        } else if (diff > 0.0) {
            s += w;
        }
    } while (w > eps || diff == 0.0);

    return val;
}


double Pull::f2(const double x, const double y) 
{
    double sump, val, tot;

    meanTmp = x;
    sigmaTmp = y;
    
    //tot = 1.0; sump = 0.2; // test
    tot = integral("fdelta", -(x_up - x_low), 0.0);
    sump = integral("fdelta", 0.0, x_up - x_low);
    sump /= (sump + tot);

    if (sump > 0.5) sump = 1.0 - sump;
    val = finverfc(sump);

    if (val >= 5.9) val = 6.01;
    return val;
}


void Pull::makeCompatPlot()
{
    if (x_low == 0.0 && x_up == 0.0 && y_low == 0.0 && y_up == 0.0) {
        x_low = copyHist->GetXaxis()->GetXmin();
        x_up = copyHist->GetXaxis()->GetXmax();
        y_low = 0.0;
        double rms = copyHist->GetRMS();
        y_up = rms * 3.0;
    }

    if (CompatPlot!=NULL) delete CompatPlot;
    CompatPlot = new TH2D("CompatPlot", "", nx, x_low, x_up, ny, y_low, y_up);

    double stepx = (x_up - x_low) / nx;
    double stepy = (y_up - y_low) / ny;
    double x, y;

    for (int i = 0; i < nx; i++) {
        x = (i + 0.5) * stepx + x_low;
        for (int j = 0; j < ny; j++) {
            y = (j + 0.5) * stepy + y_low;
            CompatPlot->Fill(x, y, f2(x, y));
        }
    }    
}



