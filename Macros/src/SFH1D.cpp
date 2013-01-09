/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <vector>
#include <numeric>
#include <TString.h>
#include <TGaxis.h>
#include <TPad.h>
#include <BAT/BCH1D.h>
#include "SFH1D.h"

SFH1D::SFH1D(TH1D& hist, const double prob68_in, const double prob95_in) 
: prob68(prob68_in), prob95(prob95_in), origHist(hist), origName(hist.GetName())
{
    std::string NewName = origName + "_new";
    newHist = (TH1D*) hist.Clone(NewName.c_str()); 
    double sum = newHist->Integral();
    newHist->Scale(1.0/sum);

    newHist68 = (TH1D*) newHist->Clone("68"); 
    newHist95 = (TH1D*) newHist->Clone("95");     
    
    computeIntervals();
}


void SFH1D::smoothHist(const int smooth)
{
    if (smooth!=0) { 
        newHist->Smooth(smooth);
        double sum = newHist->Integral();
        newHist->Scale(1.0/sum);
        computeIntervals();
    }
}


void SFH1D::increaseNbins(const int newNbin)
{
    int nbin = newHist->GetNbinsX();
    double xmin = newHist->GetXaxis()->GetXmin();
    double xmax = newHist->GetXaxis()->GetXmax();
    
    if (newNbin > nbin) {
        // compute new weights for a histogram having more bins
        std::vector<double> NewEntries;
        for (int i = 0; i < newNbin; i++) {    
            double x = xmin + (xmax - xmin)/(double)newNbin*((double)i + 0.5);
            
            int index = newHist->FindBin(x);
            if (index==1) index += 1; // the first bin
            if (index==nbin) index -= 1; // the last bin
            
            double xx[3], yy[3];
            for (int j=0; j<3; j++){
                xx[j] = newHist->GetBinCenter(index + j);
                yy[j] = newHist->GetBinContent(index + j);
            }

            double val = quadra(x, xx, yy);
            if (val < 0.0) val = 0.0;
            NewEntries.push_back(val);        
        }
        
        double sum = std::accumulate(NewEntries.begin(), NewEntries.end(), 0.0);
        
        newHist->Reset();
        newHist->SetBins(newNbin, xmin, xmax);
        for (int n = 1; n <= newNbin; n++)
            newHist->SetBinContent(n, NewEntries.at(n-1)/sum);

        computeIntervals();
    }
}


void SFH1D::Draw(const TString xlab, const TString ylab, 
                 const int col68, const int col95, const int maxDigits, 
                 const bool bOrigHist, const bool superImpose)
{
    newHist->SetTitle("");
    newHist->GetXaxis()->SetTitleSize(0.06);
    newHist->GetYaxis()->SetTitleSize(0.06);
    newHist->GetXaxis()->SetTitleOffset(1.1);
    newHist->GetYaxis()->SetTitleOffset(1.5);
    newHist->GetXaxis()->SetNdivisions(505);
    newHist->GetYaxis()->SetNdivisions(505);
    newHist->GetXaxis()->SetLabelOffset(0.013);
    newHist->GetXaxis()->SetLabelSize(0.043);
    newHist->GetYaxis()->SetLabelOffset(0.013);
    newHist->GetYaxis()->SetLabelSize(0.043);
    newHist->SetLabelFont(42,"X");
    newHist->SetLabelFont(42,"Y");
    newHist->SetTitleFont(42,"X");
    newHist->SetTitleFont(42,"Y");
    //newHist->SetLabelFont(62,"X");
    //newHist->SetLabelFont(62,"Y");
    //newHist->SetTitleFont(62,"X");
    //newHist->SetTitleFont(62,"Y");
    newHist->SetMinimum(0.0);
    newHist->SetLineStyle(1);
    ((TGaxis*) newHist->GetXaxis())->SetMaxDigits(maxDigits);
    ((TGaxis*) newHist->GetYaxis())->SetMaxDigits(maxDigits);

    // Titles of the axes 
    if (xlab.CompareTo("")==0) {
        TString Xtitle = ConvertTitle(newHist->GetXaxis()->GetTitle());
        newHist->GetXaxis()->SetTitle(Xtitle);
    } else
        newHist->GetXaxis()->SetTitle(xlab);
    if (ylab.CompareTo("")==0)
        newHist->GetYaxis()->SetTitle("Probability density");
    else 
        newHist->GetYaxis()->SetTitle(ylab);
    
    // This is necessary to set the same styles and titles as newHist 
    // to newHist68 and newHist95. 
    computeIntervals();

    // 68% and 95% intervals
    newHist68->SetFillColor(col68);
    newHist95->SetFillColor(col95);
    newHist68->SetLineStyle(3);
    newHist95->SetLineStyle(3);
    
    // normalize the histograms
    newHist->Scale(1.0/newHist->GetBinWidth(1));
    newHist68->Scale(1.0/newHist68->GetBinWidth(1));
    newHist95->Scale(1.0/newHist95->GetBinWidth(1));
            
    // draw the histograms
    if (!superImpose) 
        newHist->Draw();
    else {
        newHist->SetLineStyle(3);
        newHist68->SetFillStyle(3001);
        newHist95->SetFillStyle(3001);
    }
    newHist95->Draw("SAME");
    newHist68->Draw("SAME");
    newHist->Draw("SAME");
    
    // draw the original histogram
    if (bOrigHist) {
        origHist.SetLineColor(kRed);
        origHist.SetLineWidth(2);
        origHist.Draw("SAME");
    }

    gPad->RedrawAxis();
}


//////////////////////////////////////////////////////////////////////// 

double SFH1D::quadra(const double x, const double xx[3], const double yy[3]) const 
{
    double a = ( (yy[2]-yy[0])/(xx[2]-xx[0])*(xx[1]-xx[0]) - (yy[1]-yy[0]) )
               /( (xx[2]*xx[2]-xx[0]*xx[0])/(xx[2]-xx[0])*(xx[1]-xx[0]) 
                  - (xx[1]*xx[1]-xx[0]*xx[0]) );
    double b = ( yy[1] - yy[0] - a*(xx[1]*xx[1]-xx[0]*xx[0]) )/(xx[1]-xx[0]);
    double c = yy[0] - b*xx[0] - a*xx[0]*xx[0];
    return a*x*x + b*x + c;
}


TH1D* SFH1D::HistInterval(double &min, double &max, const double level) const
{
    BCH1D *myBCH1D = new BCH1D(newHist);  
    TH1D* tmpHist = myBCH1D->GetSmallestIntervalHistogram(level);
    for (int n = 1; n <= tmpHist->GetNbinsX(); n++) {
        if (tmpHist->GetBinContent(n)!=0.0) 
            tmpHist->SetBinContent(n, newHist->GetBinContent(n));
    }    
    
    // get the minimum of the interval
    bool flag = false;
    for (int n = 1; n <= tmpHist->GetNbinsX(); n++) {
        if (tmpHist->GetBinContent(n)!=0.0 && !flag) {
            min = tmpHist->GetBinLowEdge(n);
            flag = true;
        }
    }

    // get the maximum of the interval
    flag = false;
    for (int n = tmpHist->GetNbinsX(); n >= 1; n--) {
        if (tmpHist->GetBinContent(n)!=0.0 && !flag) {
            max = tmpHist->GetBinLowEdge(n) + tmpHist->GetBinWidth(1);
            flag = true;
        }
    }
    
    return tmpHist;
}


void SFH1D::computeIntervals()
{
    double min, max;
    if(newHist68!=NULL) delete newHist68;
    newHist68 = HistInterval(min, max, prob68);
    xmin68 = min;
    xmax68 = max;
    if(newHist95!=NULL) delete newHist95;
    newHist95 = HistInterval(min, max, prob95);
    xmin95 = min;
    xmax95 = max;
    
    localMode = newHist->GetBinCenter(newHist->GetMaximumBin());
}





