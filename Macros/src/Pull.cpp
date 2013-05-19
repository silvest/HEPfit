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
#include <Math/Functor.h>
#include <Math/IntegratorMultiDim.h>
#include <Math/AllIntegrationTypes.h>
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
    
    // draw the measured value
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


double Pull::integrand(const double* x) const
{
    double x1 = x[0];
    double y = x[1];

    double f1 = fhisto(x1);
    double f2 = TMath::Gaus(x1 - y, meanTmp, sigmaTmp, true);
    if (f1 <= 0.0) f1 = 0.0;
    if (f2 <= 0.0) f2 = 0.0;

    return ( f1*f2 );
}


double Pull::calcPull(const double mean, const double sigma, const bool lowStat)
{
    /* Mean and sigma for the indirect measurement */
    meanTmp = mean;
    sigmaTmp = sigma;

    /* parameters for numerical integrations */
    double AbsTolerance = 1.E-14; // desired absolute error (irrelevant to VEGAS)
    double RelTolerance = 1.E-6; // desired relative error (irrelevant to VEGAS)
    double ncall = 5000000; // maximum number of function calls
    if (lowStat)
        ncall = 50000; // maximum number of function calls

    /* Set the integrand which is the convolution of the two p.d.f.'s */
    ROOT::Math::Functor wf(this, &Pull::integrand, 2);

    /* The range of the histogram for the indirect measurement */
    double x1Min = copyHist->GetXaxis()->GetXmin();
    double x1Max = copyHist->GetXaxis()->GetXmax();

    /* The range of the histogram for y=x1-x2, where x1 and x2 denote the
     * random variables associated with the indirect and direct measurements
     * under consideration, respectively. */
    double yMin, yMax;
    yMin = x1Min - (meanTmp + 5.0*sigmaTmp);
    yMax = x1Max - (meanTmp - 5.0*sigmaTmp);

    /* Note that the result of the integration over the full range of x1
     is not normalized to unity, since it depends on the bin size of the
     histogram. */

    /* Total */
    ROOT::Math::IntegratorMultiDim ig2(wf, ROOT::Math::IntegrationMultiDim::kVEGAS,
                                       AbsTolerance, RelTolerance, ncall);
    double min2[2] = {x1Min, yMin};
    double max2[2] = {x1Max, yMax};
    double total = ig2.Integral(min2, max2);
    if (ig2.Status())
        std::cout << "Pull::calcPull(): Error in ig2.Integral()" << std::endl;
    if (ig2.Error() > fabs(ig2.Result())*0.01) {
        std::cout << "Pull::calcPull(): ig2.Error() > fabs(ig2.Result())*0.01" << std::endl;
        std::cout << "total= " << total << " +- " << ig2.Error()
                  << " [status=" << ig2.Status() << "]" << std::endl;
    }

    ROOT::Math::IntegratorMultiDim ig(wf, ROOT::Math::IntegrationMultiDim::kVEGAS,
                                      AbsTolerance, RelTolerance, ncall);
    double val;
    if (0.0 > yMax)
        val = 0.0;
    else if (0.0 < yMin)
        val = 1.0;
    else {
        if ( copyHist->GetMean() - meanTmp < 0.0 ) {
            double min[2] = {x1Min, 0.0};
            double max[2] = {x1Max, yMax};
            val = ig.Integral(min, max)/total;
        } else {
            double min[2] = {x1Min, yMin};
            double max[2] = {x1Max, 0.0};
            val = 1.0 - ig.Integral(min, max)/total;
        }
        if (ig.Status())
            std::cout << "Pull::calcPull(): Error in ig.Integral()" << std::endl;
        if (ig.Error() > fabs(ig.Result())*0.01) {
            std::cout << "Pull::calcPull(): ig.Error() > fabs(ig.Result())*0.01" << std::endl;
            std::cout << "val= " << ig.Result() << " +- " << ig.Error()
                      << " [status=" << ig.Status() << "]" << std::endl;
        }
    }


    /* Pull in units of the standard deviation */
    double pull = 0.0;
    if ( fabs(1.0-2.0*val) < 1.0  ) {
        /* ErfInverse(x): x must be -1<x<1 */
        //pull = - sqrt(2.0)*TMath::ErfInverse(1.0 - 2.0*val);
        /* ErfcInverse(x): x must be 0<x<2 */
        pull = - sqrt(2.0)*TMath::ErfcInverse(2.0*val);
    } else {
        double sign = 1.0;
        if (1.0-2.0*val < 0.0) sign = -1.0;
        pull = - 6.01 * sign;
    }

    return pull;
}


void Pull::makeCompatPlot()
{
    if (x_low == 0.0 && x_up == 0.0 && y_low == 0.0 && y_up == 0.0) {
        x_low = copyHist->GetXaxis()->GetXmin();
        x_up = copyHist->GetXaxis()->GetXmax();
        y_low = 0.0;
        y_up = copyHist->GetRMS() * 3.0;
    }

    if (CompatPlot!=NULL) delete CompatPlot;
    CompatPlot = new TH2D("CompatPlot", "", nx, x_low, x_up, ny, y_low, y_up);

    double stepx = (x_up - x_low) / (double)nx;
    double stepy = (y_up - y_low) / (double)ny;
    double x, y;

    double pull;
    std::cout << "nx: " << std::endl;
    for (int i = 0; i < nx; i++) {
        if (i%10==0.0) std::cout << i << "/" << nx << std::endl;
        x = ((double)i + 0.5) * stepx + x_low;
        for (int j = 0; j < ny; j++) {
            y = ((double)j + 0.5) * stepy + y_low;
            pull = fabs(calcPull(x, y, true));
            if (pull > 6.0) pull = 6.01;

            CompatPlot->Fill(x, y, pull);
        }
    }
}



