/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <vector>
#include <algorithm>
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
    newHist->Scale(1.0/newHist->Integral());

    HistAxes = NULL;
    
    newHist68 = (TH1D*) newHist->Clone("68"); 
    newHist95 = (TH1D*) newHist->Clone("95");     

    computeIntervals();
}


void SFH1D::smoothHist(const int smooth)
{
    if (smooth!=0) { 
        newHist->Smooth(smooth);
        newHist->Scale(1.0/newHist->Integral());
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


void SFH1D::DrawAxes(const TString xlab, const TString ylab, 
                     const int maxDigits, 
                     const double x_low, const double x_up)
{    
    // the range of the x axis    
    double xLow, xUp;
    if (x_low == 0.0 && x_up == 0.0) {
        xLow = newHist->GetXaxis()->GetXmin();
        xUp = newHist->GetXaxis()->GetXmax();
    } else {
        xLow = x_low;
        xUp = x_up;
    }

    newHist->Scale(1.0/newHist->Integral());
    HistAxes = (TH1D*) newHist->Clone();    
    HistAxes->Reset("M");
    delete HistAxes;
    
    HistAxes = new TH1D("HistAxes","HistAxes", newHist->GetNbinsX(), xLow, xUp); 

    // the range of the y axis
    HistAxes->SetMinimum(0.0);
    HistAxes->SetMaximum(1.1*newHist->GetMaximum()/newHist->GetBinWidth(1));
    
    HistAxes->SetTitle("");
    HistAxes->GetXaxis()->SetTitleSize(0.06);
    HistAxes->GetYaxis()->SetTitleSize(0.06);
    HistAxes->GetXaxis()->SetTitleOffset(1.1);
    HistAxes->GetYaxis()->SetTitleOffset(1.5);
    HistAxes->GetXaxis()->SetNdivisions(505);
    HistAxes->GetYaxis()->SetNdivisions(505);
    HistAxes->GetXaxis()->SetLabelOffset(0.013);
    HistAxes->GetXaxis()->SetLabelSize(0.043);
    HistAxes->GetYaxis()->SetLabelOffset(0.013);
    HistAxes->GetYaxis()->SetLabelSize(0.043);
    HistAxes->SetLabelFont(42,"X");
    HistAxes->SetLabelFont(42,"Y");
    HistAxes->SetTitleFont(42,"X");
    HistAxes->SetTitleFont(42,"Y");
    //HistAxes->SetLabelFont(62,"X"); // bold font
    //HistAxes->SetLabelFont(62,"Y");
    //HistAxes->SetTitleFont(62,"X");
    //HistAxes->SetTitleFont(62,"Y");
    ((TGaxis*) HistAxes->GetXaxis())->SetMaxDigits(maxDigits);
    ((TGaxis*) HistAxes->GetYaxis())->SetMaxDigits(maxDigits);
    
    // Titles of the axes 
    if (xlab.CompareTo("")==0) {
        TString Xtitle = ConvertTitle(newHist->GetXaxis()->GetTitle());
        HistAxes->GetXaxis()->SetTitle(Xtitle);
    } else
        HistAxes->GetXaxis()->SetTitle(xlab);
    if (ylab.CompareTo("")==0)
        HistAxes->GetYaxis()->SetTitle("Probability density");
    else 
        HistAxes->GetYaxis()->SetTitle(ylab);
    
    HistAxes->Draw();    
}


void SFH1D::Draw(const int lineStyle, const int lineWidth, const int lineColor, 
                 const int col68, const int col95, const int fillStyle, 
                 const bool bOrigHist)
{    
    // normalize the histograms
    newHist->Scale(1.0/newHist->Integral());
    newHist68->Scale(1.0/newHist->Integral());
    newHist95->Scale(1.0/newHist->Integral());
    newHist->Scale(1.0/newHist->GetBinWidth(1));
    newHist68->Scale(1.0/newHist->GetBinWidth(1));
    newHist95->Scale(1.0/newHist->GetBinWidth(1));
    
    // styles
    newHist->SetLineStyle(lineStyle);
    newHist68->SetLineStyle(3);
    newHist95->SetLineStyle(3);
    newHist->SetLineWidth(lineWidth);
    newHist68->SetLineWidth(lineWidth);
    newHist95->SetLineWidth(lineWidth);
    newHist->SetLineColor(lineColor);
    newHist68->SetLineColor(lineColor);
    newHist95->SetLineColor(lineColor);
    newHist68->SetFillColor(col68);
    newHist95->SetFillColor(col95);
    newHist68->SetFillStyle(fillStyle);
    newHist95->SetFillStyle(fillStyle);

    if (prob95 != 0.0) newHist95->Draw("SAME");
    if (prob68 != 0.0) newHist68->Draw("SAME");
    newHist->Draw("SAME");
    
    // draw the original histogram
    if (bOrigHist) {
        origHist.SetLineColor(kRed);
        origHist.SetLineWidth(2);
        origHist.Draw("SAME");
    }

    gPad->RedrawAxis();
}


void SFH1D::RescaleYaxis(double max_another)
{
    double max_new = std::max( HistAxes->GetMaximum(), 1.1*max_another );
    HistAxes->SetMaximum(max_new);    
}


void SFH1D::OutputResults(ostream& os, const int smooth, const bool WasDrawed) const 
{
    //   Note: after Draw(), use Integral("width").     
    char opt[10] = "";
    if(WasDrawed) strcpy(opt, "width");
    
    os << "  Num of bins: " << newHist->GetNbinsX() 
       << "   smooth: " << smooth << " time(s)" << std::endl
       << "  Local mode: " << localMode
       << " + " << xmax68 - localMode 
       << " - " << localMode - xmin68 << std::endl
       << "  Center of " << newHist68->Integral(opt)*100.0 
       << "% interval: " << (xmin68 + xmax68)/2.0
       << " +- " << (xmax68 - xmin68)/2.0 << std::endl
       << "  at " << newHist68->Integral(opt)*100.0 << " % (>~"
       << prob68*100.0 << "%)" << " [" << xmin68 << ", " << xmax68 << "]" 
       << std::endl
       << "  at " << newHist95->Integral(opt)*100.0 << " % (>~"
       << prob95*100.0 << "%)" << " [" << xmin95 << ", " << xmax95 << "]" 
       << std::endl;        
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
    newHist->Scale(1.0/newHist->Integral());
    
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





