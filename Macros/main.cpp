/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>
#include <TObject.h>
#include <TF1.h>
#include <TLegend.h>
#include "BaseMacros.h"
#include "SFH1D.h"
#include "Pull.h"
#include "SFH2D.h"

using namespace std;

/**
 *  @addtogroup Macros
 *  A project for drawing histograms and outputting results. 
 *  @{
 */

//-- set probability ranges --
//const double prob68 = 0.6827;
//const double prob95 = 0.9545;
const double prob68 = 0.68;
const double prob95 = 0.95;

/**
 * @brief 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
int main(int argc, char** argv) 
{
    
    if (argc < 3) {
        cout << "#################################################################" << endl;
        cout << " Usage:                                                          " << endl;
        cout << "   1-D histogram:                                                " << endl;
        cout << "     > " << argv[0] << " rootfile plotname --oneDim <optional parameters>" << endl;
        cout << "   Compatibility plot:                                           " << endl;
        cout << "     > " << argv[0] << " rootfile plotname --compat <optional parameters>" << endl;
        cout << "   2-D histogram:                                                " << endl;
        cout << "     > " << argv[0] << " rootfile plotname --twoDim <optional parameters>" << endl;
        cout << "       'rootfile' (with extension)                               " << endl;
        cout << "       'plotname' (e.g. sin2b)                                   " << endl;
        cout << "                                                                 " << endl;
        cout << "   Check the contents of a rootfile:                             " << endl;
        cout << "     > " << argv[0] << " rootfile info                           " << endl;
        cout << "                                                                 " << endl;
        cout << " Optional parameters for 1-D histograms:                         " << endl;
        cout << "   --oneDim         -> 1-D histogram (mandatory)                 " << endl;
        cout << "   --orig           -> superimpose the original histogram        " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   -smooth=ntime    -> iterative smoothing with TH1::Smooth()    " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -moreBins=nbin      -> increase the number of bins            " << endl;
        cout << "                       (default: nbin=100)                       " << endl;
        cout << "   -col68=index     -> color index of the 68% interval           " << endl;
        cout << "                       (default: index=1393)                     " << endl;
        cout << "   -col95=index     -> color index of the 95% interval           " << endl;
        cout << "                       (default: index=1392)                     " << endl;
        cout << "   -plot2=name      -> plot to be superimposed                   " << endl;
        cout << "   -rootfile2=filename2 -> rootfile name for plot2 (with extension)" << endl;
        cout << "                           (default: same as rootfile)           " << endl;
        cout << "   -smooth2=ntime   -> iterative smoothing for plot2             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -moreBins2=nbin  -> increase the number of bins for plot2     " << endl;
        cout << "                       (default: nbin=100)                       " << endl;
        cout << "   -col682=index    -> color index of the 68% interval for plot2 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col952=index    -> color index of the 95% interval for plot2 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -priorMean=mean  -> mean value for a Gaussian function to be superimposed" << endl;
        cout << "   -priorSigma=sig  -> standard deviation for a Gaussian function to be superimposed" << endl;
        cout << "                                                                 " << endl;
        cout << " Optional parameters for compatibility plots:                    " << endl;
        cout << "   --compat           -> Compatibility plot (mandatory)          " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   -val=xval        -> xval is the measured point to be superimposed" << endl;
        cout << "   -err=xerr        -> xerr is the error of xval                 " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range  " << endl;
        cout << "   -bins=\"[xbins]x[ybins]\"          -> define the graph binning" << endl;
        cout << "                                         (default: 100x20)       " << endl;
        cout << "                                                                 " << endl;
        cout << " Optional parameters for 2-D histograms:                         " << endl;
        cout << "   --twoDim         -> 2-D histogram (mandatory)                 " << endl;
        cout << "   --drawlines      -> draw contour lines                        " << endl;
        cout << "   --outputTxt      -> output results to a text file             " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       (default: plotname)                       " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels         " << endl;
        cout << "                       (default: n=8)                            " << endl;
        cout << "   -precision=n     -> precision of values in the std output     " << endl;
        cout << "                       (default: n=6)                            " << endl;
        cout << "   -xlab=namex      -> x label                                   " << endl;
        cout << "   -ylab=namey      -> y label                                   " << endl;
        cout << "   -smooth=ntime    -> iterative smoothing with TH2::Smooth()    " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col68=index     -> color index of the 68% interval           " << endl;
        cout << "                       (default: index=1393)                     " << endl;
        cout << "   -col95=index     -> color index of the 95% interval           " << endl;
        cout << "                       (default: index=1392)                     " << endl;
        cout << "   -plot2=name      -> plot to be superimposed                   " << endl;
        cout << "   -rootfile2=filename -> rootfile name for plot2 (with extension)" << endl;
        cout << "                          (default: same as rootfile)            " << endl;
        cout << "   -smooth2=ntime   -> iterative smoothing for plot2             " << endl;
        cout << "                       (default: ntime=0)                        " << endl;
        cout << "   -col682=index    -> color index of the 68% interval for plot2 " << endl;
        cout << "                       (default: index=2)                        " << endl;
        cout << "   -col952=index    -> color index of the 95% interval for plot2 " << endl;
        cout << "                       (default: index=5)                        " << endl;
        cout << "   -Xval=xval2      -> xval2 is the  point to be superimposed    " << endl;
        cout << "   -Xerr=xerr2      -> xerr2 is the error of xval2               " << endl;
        cout << "   -Yval=yval2      -> yval2 is the  point to be superimposed    " << endl;
        cout << "   -Yerr=yerr2      -> yerr2 is the error of yval2               " << endl;
        cout << "#################################################################" << endl;
        return 0;
    }

    // parameters which can be changed by command-line arguments
    bool bOneDim = false, bCompat = false, bTwoDim = false;
    bool bOrig = false, bOutputTxt = false, bMoreBins = false;
    bool bContLines = false, bSuperImpose = false;
    int maxDig = 8, prec = 6, smooth = 0, newNbins = 100, col68 = 1393, col95 = 1392;
    int nx = 100, ny = 20, smooth2 = 0, newNbins2 = 100, col682 = 2, col952 = 5;
    double xval = -999.0, xerr = 0.0, x_low = 0.0, x_up = 0.0, y_low = 0.0, y_up = 0.0;
    double xval2 = -999.0, xerr2 = 0.0, yval2 = -999.0, yerr2 = 0.0;
    double prior_mean = 0.0, prior_sigma = 0.0;
    string plotname2, filename2;
    TString xlab = "", ylab = "";

    // root file
    string filename = argv[1];
    TFile *datafile = new TFile(filename.c_str());
    if (argc == 3 && strncmp(argv[2], "info", 4) == 0) {
        datafile->ls();
        return 0;
    }
        
    // plot name
    string plotname = argv[2];
    TObject* tobj;
    tobj = datafile->Get(plotname.c_str());
    if (tobj == NULL) {
        cout << "Error: plot \"" << plotname << "\" does not exist..." << endl;
        return 1;
    }
    
    // options
    string outputFileName = plotname;
    filename2 = filename;
    for (int i = 3; i < argc; i++) {
        char str[100];
        if (strncmp(argv[i], "--oneDim", 8) == 0) bOneDim = true;
        if (strncmp(argv[i], "--compat", 8) == 0) bCompat = true;
        if (strncmp(argv[i], "--twoDim", 8) == 0) bTwoDim = true;    
        if (strncmp(argv[i], "--orig", 6) == 0) bOrig = true;
        if (strncmp(argv[i], "--outputTxt", 11) == 0) bOutputTxt = true;        
        if (strncmp(argv[i], "--drawlines", 11) == 0) bContLines = true;        
        
        if (strncmp(argv[i], "-output", 7) == 0) {
            sscanf(argv[i], "-output=%s", str);
            outputFileName = str;
        }

        if (strncmp(argv[i], "-xlab", 5) == 0) {
            sscanf(argv[i], "-xlab=%s", str);
            xlab = str;
        }
        
        if (strncmp(argv[i], "-ylab", 5) == 0) {
            sscanf(argv[i], "-ylab=%s", str);
            ylab = str;
        }        
        
        if (strncmp(argv[i], "-maxDigits", 10) == 0) 
            sscanf(argv[i], "-maxDigits=%d", &maxDig);

        if (strncmp(argv[i], "-precision", 10) == 0) 
            sscanf(argv[i], "-precision=%d", &prec);

        if (strncmp(argv[i], "-smooth", 7) == 0)
            sscanf(argv[i], "-smooth=%d", &smooth);
            
        if (strncmp(argv[i], "-smooth2", 8) == 0) 
            sscanf(argv[i], "-smooth2=%d", &smooth2);

        if (strncmp(argv[i], "-moreBins", 9) == 0) {
            sscanf(argv[i], "-moreBins=%d", &newNbins);
            bMoreBins = true;
        }
        
        if (strncmp(argv[i], "-moreBins2", 10) == 0) {
            sscanf(argv[i], "-moreBins2=%d", &newNbins2);
            bMoreBins = true;
        }

        if (strncmp(argv[i], "-col68", 6) == 0) 
            sscanf(argv[i], "-col68=%d", &col68);
        
        if (strncmp(argv[i], "-col95", 6) == 0) 
            sscanf(argv[i], "-col95=%d", &col95);
        
        if (strncmp(argv[i], "-priorMean", 10) == 0) 
            sscanf(argv[i], "-priorMean=%lf", &prior_mean);
        
        if (strncmp(argv[i], "-priorSigma", 11) == 0) 
            sscanf(argv[i], "-priorSigma=%lf", &prior_sigma);

        if (strncmp(argv[i], "-val", 4) == 0) 
            sscanf(argv[i], "-val=%lf", &xval);

        if (strncmp(argv[i], "-err", 4) == 0) 
            sscanf(argv[i], "-err=%lf", &xerr);

        if (strncmp(argv[i], "-range", 6) == 0) {
            TString stmp(argv[i] + 7);
            sscanf(stmp.Data(), "[%lf,%lf]x[%lf,%lf]", &x_low, &x_up, &y_low, &y_up);
        }

        if (strncmp(argv[i], "-bins", 5) == 0) {
            TString stmp(argv[i] + 6);
            sscanf(stmp.Data(), "[%d]x[%d]", &nx, &ny);
        }
        
        // second plot
        if (strncmp(argv[i], "-plot2", 6) == 0) {
            sscanf(argv[i], "-plot2=%s", str);
            plotname2 = str;
            bSuperImpose = true;
        }

        // root file for the second plot
        if (strncmp(argv[i], "-rootfile2", 10) == 0) {
            sscanf(argv[i], "-rootfile2=%s", str);
            filename2 = str;
        }        
        
        if (strncmp(argv[i], "-col682", 7) == 0) 
            sscanf(argv[i], "-col682=%d", &col682);
        
        if (strncmp(argv[i], "-col952", 7) == 0) 
            sscanf(argv[i], "-col952=%d", &col952);
        
        if (strncmp(argv[i], "-Xval", 5) == 0) 
            sscanf(argv[i], "-Xval=%lf", &xval2);

        if (strncmp(argv[i], "-Xerr", 5) == 0) 
            sscanf(argv[i], "-Xerr=%lf", &xerr2);

        if (strncmp(argv[i], "-Yval", 5) == 0) 
            sscanf(argv[i], "-Yval=%lf", &yval2);

        if (strncmp(argv[i], "-Yerr", 5) == 0) 
            sscanf(argv[i], "-Yerr=%lf", &yerr2);
    }        

    if (!(bOneDim | bCompat | bTwoDim)) {
        cout << "Error: A mandatory option --oneDim, --compat or --twoDim is missing!" << endl;
        return 1;
    }
    if ( (bOneDim&bCompat) | (bOneDim&bTwoDim) | (bTwoDim&bCompat) ) {
        cout << "Error: Use only one of --oneDim, --compat and --twoDim!" << endl;
        return 1;
    }

    string epsFileName = outputFileName + ".eps";
    string txtFileName = outputFileName + ".txt";

    // second plot
    TFile *datafile2 = new TFile(filename2.c_str());
    TObject* tobj2;
    tobj2 = datafile2->Get(plotname2.c_str());
    if (tobj2 == NULL && bSuperImpose) {
        cout << "Error: plot \"" << plotname2 << "\" does not exist..." << endl;
        return 1;
    }
    
    //----------------------------------------------------------------------

    ofstream* fout;
    streambuf* buf;
    if (bOutputTxt) {
        fout = new ofstream();
        fout->open(txtFileName.c_str(), ios::out);
        if (!fout->is_open())
            cout << "cannot open " << txtFileName << endl;
        buf = fout->rdbuf();        
    } else
        buf = cout.rdbuf();
    ostream os(buf);
    os.precision(prec);

    BaseMacros::DefineNewColours();
    
    TCanvas TC("TC", "2D histogram", 3);
    TC.SetLeftMargin(0.20);
    // TC.SetRightMargin(0.20);
    TC.SetBottomMargin(0.15);
    TC.SetFillColor(0);
    TC.SetBorderMode(0);
    TC.SetFrameFillColor(0);
    TC.SetHighLightColor(0);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    //gStyle->SetStripDecimals(false);   
                
    if (bOneDim) {
        // 1-D histogram
        TH1D* hist = (TH1D*) tobj->Clone();   
        os << hist->GetXaxis()->GetTitle() << " in " << plotname << endl;

        // output the 1-D histogram
        SFH1D SFHisto1D(*hist, prob68, prob95);
        SFHisto1D.smoothHist(smooth);
        SFHisto1D.increaseNbins(newNbins);
        SFHisto1D.Draw(xlab, ylab, col68, col95, maxDig, bOrig, false);    
    
        // rescale
        //SFHisto1D.getNewHist()->Scale(10.0);
        //SFHisto1D.getNewHist68()->Scale(10.0);
        //SFHisto1D.getNewHist95()->Scale(10.0);
        //SFHisto1D.getNewHist()->GetXaxis()->SetRange(400,1400);
        
        // superimpose another 1-D histogram (e.g. for a posterior)
        // and a Gaussian (prior) function
        if (prior_sigma != 0.0) {
            // another 1-D histogram
            double ymax_new;
            TH1D* hist2 = (TH1D*) tobj2->Clone();
            SFH1D SFHisto1D2(*hist2, prob68, prob95);
            if (bSuperImpose) {
                SFHisto1D2.smoothHist(smooth2);
                SFHisto1D2.increaseNbins(newNbins2);
                SFHisto1D2.Draw("", "", col682, col952, maxDig, bOrig, bSuperImpose);    
                ymax_new = max( SFHisto1D.getNewHist()->GetMaximum(), 
                                1.1*SFHisto1D2.getNewHist()->GetMaximum() );
                SFHisto1D.getNewHist()->SetMaximum(ymax_new);
            }
            
            // a Gaussian function 
            TF1* prior = new TF1("prior",
                    "1./sqrt(2.*TMath::Pi())/[1]* exp(- (x-[0])*(x-[0])/2./[1]/[1])",
                    SFHisto1D.getNewHist()->GetXaxis()->GetXmin(),
                    SFHisto1D.getNewHist()->GetXaxis()->GetXmax());
            prior->SetParameter(0, prior_mean);
            prior->SetParameter(1, prior_sigma);    
            prior->SetLineStyle(2);
            prior->SetLineWidth(4);
            prior->SetNpx(1000);
            ymax_new = max( 1.1*prior->GetMaximum(), 
                            SFHisto1D.getNewHist()->GetMaximum() );
            SFHisto1D.getNewHist()->SetMaximum(ymax_new);
            prior->Draw("SAME");

            // draw the legend
            TLegend *legend = new TLegend(0.63,0.65,0.88,0.85);
            //TLegend *legend = new TLegend(0.23,0.65,0.48,0.85);
            legend->SetFillColor(0);
            legend->SetBorderSize(0);
            legend->SetTextFont(42);
            legend->SetTextSize(0.043);
            legend->SetMargin(0.4);
            legend->AddEntry(prior,"Prior","L");
            if (bSuperImpose) 
                legend->AddEntry(SFHisto1D2.getNewHist68(),"Posterior","F");
            legend->AddEntry(SFHisto1D.getNewHist68(),"Fit","F");
            //legend->AddEntry(SFHisto1D.getNewHist68(),"Fit [x10]","F");
            legend->Draw("");
        }
        
        // output results
        //   Note: after Draw(), use Integral("width"). 
        os << "  Num of bins: " << SFHisto1D.getNewHist()->GetNbinsX() 
           << "   smooth: " << smooth << " time(s)" << endl
           << "  Local mode: " << SFHisto1D.getLocalMode()
           << " + " << SFHisto1D.getXmax68() - SFHisto1D.getLocalMode() 
           << " - " << SFHisto1D.getLocalMode() - SFHisto1D.getXmin68() << endl
           << "  Center of " << SFHisto1D.getNewHist68()->Integral("width")*100.0 
           << "% interval: " << (SFHisto1D.getXmin68() + SFHisto1D.getXmax68())/2.0
           << " +- " << (SFHisto1D.getXmax68() - SFHisto1D.getXmin68())/2.0 << endl
           << "  at " << SFHisto1D.getNewHist68()->Integral("width")*100.0 << " % (>~"
           << prob68*100.0 << "%)"
           << " [" << SFHisto1D.getXmin68() << ", " << SFHisto1D.getXmax68() << "]" 
           << endl
           << "  at " << SFHisto1D.getNewHist95()->Integral("width")*100.0 << " % (>~"
           << prob95*100.0 << "%)"
           << " [" << SFHisto1D.getXmin95() << ", " << SFHisto1D.getXmax95() << "]" 
           << endl;        
        SFH1D orig1D(*hist, prob68, prob95);
        os << endl << "[Original histogram]" << endl
           << "  Num of bins: " << orig1D.getNewHist()->GetNbinsX() 
           << "   smooth: 0 time" << endl
           << "  Local mode: " << orig1D.getLocalMode()
           << " + " << orig1D.getXmax68() - orig1D.getLocalMode() 
           << " - " << orig1D.getLocalMode() - orig1D.getXmin68() << endl
           << "  Center of " << orig1D.getNewHist68()->Integral("")*100.0 
           << "% interval: " << (orig1D.getXmin68() + orig1D.getXmax68())/2.0
           << " +- " << (orig1D.getXmax68() - orig1D.getXmin68())/2.0 << endl
           << "  at " << orig1D.getNewHist68()->Integral("")*100.0 << " % (>~"
           << prob68*100.0 << "%)"
           << " [" << orig1D.getXmin68() << ", " << orig1D.getXmax68() << "]" 
           << endl
           << "  at " << orig1D.getNewHist95()->Integral("")*100.0 << " % (>~"
           << prob95*100.0 << "%)"
           << " [" << orig1D.getXmin95() << ", " << orig1D.getXmax95() << "]" 
           << endl;          
    } else if (bCompat) {
        // 1-D histogram
        TH1D* hist = (TH1D*) tobj->Clone();
        os << hist->GetXaxis()->GetTitle() << " in " << plotname << endl;
        os << "  Num of bins: " << nx << " x " << ny << endl;
        
        // output the compatibility plot
        Pull CompatPlot(*hist, nx, ny, x_low, x_up, y_low, y_up);
        CompatPlot.Draw(xlab, ylab, xval, xerr, maxDig);
        
    } else if (bTwoDim) {    
        // 2-D histogram
        TH2D* hist = (TH2D*) tobj->Clone();
        os << hist->GetXaxis()->GetTitle() << " vs " << hist->GetYaxis()->GetTitle() 
           << " in " << plotname << endl;
        os << "  smooth: " << smooth << " time(s)" << endl;
        
        // output the 2-D histogram
        SFH2D SFHisto2D(*hist, os, prob68, prob95);
        SFHisto2D.smoothHist(smooth);
        SFHisto2D.draw(xlab, ylab, col68, col95, maxDig, 
                       xval2, xerr2, yval2, yerr2, bContLines, false);

        // superimpose the 2nd histogram
        if (bSuperImpose) {
            TH2D* hist2 = (TH2D*) tobj2->Clone();
            SFH2D SFHisto2D2(*hist2, os, prob68, prob95);
            SFHisto2D2.smoothHist(smooth2);
            SFHisto2D2.draw(xlab, ylab, col682, col952, maxDig, 
                            xval2, xerr2, yval2, yerr2, bContLines, true);
        }
    } 
    
    gPad->RedrawAxis();
    TC.Print(epsFileName.c_str());
    
    if (bOutputTxt) {
        fout->close();          
        cout << txtFileName << " has been created." << endl;
    }

    return 0;
}

/** 
 * @}
 */


