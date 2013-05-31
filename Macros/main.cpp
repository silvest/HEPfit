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
#include <sstream>
#include <algorithm>
#include <TRoot.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TColor.h>
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
 * @addtogroup Macros
 * @brief A project for drawing histograms and outputting results. 
 * @{
 */

/**
 * @brief 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
int main(int argc, char** argv) 
{

    if (argc < 3) {
        cout << "######################################################################" << endl;
        cout << " Usage:                                                               " << endl;
        cout << "   1-D histogram:                                                     " << endl;
        cout << "     > macros --oneDim -rootfile1=rootfile -plot1=plotname <optional parameters>" << endl;
        cout << "   Compatibility plot:                                                " << endl;
        cout << "     > macros --compat -rootfile1=rootfile -plot1=plotname <optional parameters>" << endl;
        cout << "   2-D histogram:                                                     " << endl;
        cout << "     > macros --twoDim -rootfile1=rootfile -plot1=plotname <optional parameters>" << endl;
        cout << "       'rootfile' (with extension)                                    " << endl;
        cout << "       'plotname' (e.g. sin2b)                                        " << endl;
        cout << "                                                                      " << endl;
        cout << "   Check the contents of a rootfile:                                  " << endl;
        cout << "     > macros info -rootfile1=rootfile                                " << endl;
        cout << "                                                                      " << endl;
        cout << " Optional parameters for 1-D histograms:                              " << endl;
        cout << "   --oneDim         -> 1-D histogram (mandatory)                      " << endl;
        cout << "   --orig           -> superimpose the original histogram             " << endl;
        cout << "   --leftLegend     -> put the legend in the left area                " << endl;
        cout << "   --outputTxt      -> output results to a text file                  " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       [default: plotname]                            " << endl;
        cout << "   --pdf            -> output pdf instead of eps                      " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -prob68=prob     -> probability for the first interval [default: prob=0.68]" << endl;
        cout << "   -prob95=prob     -> probability for the second interval [default: prob=0.95]"<< endl;
        cout << "   -prob99=prob     -> probability for the third interval [default: prob=0.99]"<< endl;
        cout << "   -xlab=namex      -> title of the x axis                            " << endl;
        cout << "   -ylab=namey      -> title of the y axis                            " << endl;
        cout << "   -yTitleOffset=d  -> offset of the title of the y axis [default: d=1.5]" << endl; 
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -xrange=\"[xmin,xmax]\" -> define the graph range                  " << endl;
        cout << "   -legScale=scale  -> scale factor for the legend [default: scale=1] " << endl;
        cout << "   -legXmin=xmin    -> xmin for the position of the legend [default: xmin=0.63]" << endl;
        cout << "   -legYmax=ymax    -> ymax for the position of the legend [default: ymax=0.88]" << endl;
        cout << "   *** options for histograms (N=1,2,3,4,5) ***                       " << endl;
        cout << "   -plotN=name      -> name of the histogram                          " << endl;
        cout << "   -rootfileN=filename -> rootfile name for plotN (with extension)    " << endl;
        cout << "         [For N=2,3, the same as the rootfile for plot1 by default]   " << endl;
        cout << "   --onlyLineN      -> draw only the line                             " << endl;
        cout << "   -smoothN=ntime   -> iterative smoothing for plotN [default: ntime=0]" << endl;
        cout << "   -moreBinsN=nbin  -> increase the number of bins for plotN [default: nbin=100]" << endl;
        cout << "   -lineStyleN=index -> index of the line style [default: index=1]    " << endl;
        cout << "   -lineWidthN=index -> index of the line width [default: index=1]    " << endl;
        cout << "   -lineColorN=index -> index of the line color [default: index=1]    " << endl;
        cout << "   -col68N=index    -> color index of the 68% interval for plotN [default: index=1393]" << endl;
        cout << "   -col95N=index    -> color index of the 95% interval for plotN [default: index=1392]" << endl;
        cout << "   -col68alphaN=alpha -> color transparency of the 68% interval for plotN [default: alpha=1]" << endl;
        cout << "   -col95alphaN=alpha -> color transparency of the 95% interval for plotN [default: alpha=1]" << endl;
        cout << "   -fillN=index     -> index of the fill area style [default: index=1001]" << endl;
        cout << "   -legN=legend     -> legend for plotN [default: no legend]          " << endl;
        cout << "   [e.g. -plot2=Mw -roogfile2=test.root -smooth2=3 -col682=2 -col952=5]" << endl;
        cout << "   *** superimpose a Gaussian function ***                            " << endl;
        cout << "   -priorMean=mean  -> mean value for the Gaussian function           " << endl;
        cout << "   -priorSigma=sig  -> standard deviation for the Gaussian function   " << endl;
        cout << "   -legGauss=legend -> legend for the Gaussian function [default: no legend]" << endl;
        cout << "   --rescaleForMHl  -> rescale the histogram for mHl                  " << endl;
        cout << "                                                                      " << endl;
        cout << " Optional parameters for compatibility plots:                         " << endl;
        cout << "   --compat         -> Compatibility plot (mandatory)                 " << endl;
        cout << "   --outputTxt      -> output results to a text file                  " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       [default: plotname]                            " << endl;
        cout << "   --pdf            -> output pdf instead of eps                      " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -prob68=prob     -> probability for the first interval [default: prob=0.68]" << endl;
        cout << "   -prob95=prob     -> probability for the second interval [default: prob=0.95]"<< endl;
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range       " << endl;
        cout << "   -bins=\"[xbins]x[ybins]\"  -> define the graph binning [default: 100x20]" << endl;
        cout << "   -xlab=namex      -> title of the x axis                            " << endl;
        cout << "   -ylab=namey      -> title of the y axis                            " << endl;
        cout << "   -yTitleOffset=d  -> offset of the title of the y axis [default: d=1.5]" << endl; 
        cout << "   *** put a measured point ***                                       " << endl;
        cout << "   -val=xval        -> xval is the measured point                     " << endl;
        cout << "   -err=xerr        -> xerr is the error of xval                      " << endl;
        cout << "                                                                      " << endl;
        cout << " Optional parameters for 2-D histograms:                              " << endl;
        cout << "   --twoDim         -> 2-D histogram (mandatory)                      " << endl;
        cout << "   --drawlines      -> draw contour lines                             " << endl;
        cout << "   --outputTxt      -> output results to a text file                  " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       [default: plotname]                            " << endl;
        cout << "   --pdf            -> output pdf instead of eps                      " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range       " << endl;
        cout << "   -xlab=namex      -> title of the x axis                            " << endl;
        cout << "   -ylab=namey      -> title of the y axis                            " << endl;
        cout << "   -yTitleOffset=d  -> offset of the title of the y axis [default: d=1.5]" << endl; 
        cout << "   -legScale=scale  -> scale factor for the legend [default: scale=1.0]" << endl;
        cout << "   -legXmin=xmin    -> xmin for the position of the legend [default: xmin=0.63]" << endl;
        cout << "   -legYmax=ymax    -> ymax for the position of the legend [default: ymax=0.88]" << endl;
        cout << "   *** options for histogram (N=1,2,3,4,5) ***                        " << endl;
        cout << "   --only95N        -> draw only the 95% contour                      " << endl;
        cout << "   -plotN=name      -> name of the histogram                          " << endl;
        cout << "   -rootfileN=filename -> rootfile name for plotN (with extension)    " << endl;
        cout << "         [For N=2,3, the same as the rootfile for plot1 by default]   " << endl;
        cout << "   -smoothN=ntime   -> iterative smoothing for plotN [default: ntime=0]" << endl;
        cout << "   -col68N=index    -> color index of the 68% interval for plotN [default: index=1393]" << endl;
        cout << "   -col95N=index    -> color index of the 95% interval for plotN [default: index=1392]" << endl;
        cout << "   -col68alphaN=alpha -> color transparency of the 68% interval for plotN [default: alpha=1]" << endl;
        cout << "   -col95alphaN=alpha -> color transparency of the 95% interval for plotN [default: alpha=1]" << endl;
        cout << "   -lineStyleN=index-> index of the line style [default: index=1]     " << endl;
        cout << "   -lineStyle68N=index-> index of the line style for the 68% contour [default: index=3]" << endl;
        cout << "   -lineWidthN=index-> index of the line width [default: index=1]    " << endl;
        cout << "   -lineColorN=index-> index of the line color [default: index=1]    " << endl;
        cout << "   -fillN=index     -> index of the fill area style [default: index=1001]" << endl;
        cout << "   -legN=legend     -> legend for plotN [default: no legend]          " << endl;
        cout << "   -NumNewPointsN=num-> add more points to close contours [default: num=20]" << endl;
        cout << "   *** put a measured point ***                                       " << endl;
        cout << "   -Xval=xval2      -> xval2 is the x value of the point              " << endl;
        cout << "   -Xerr=xerr2      -> xerr2 is the error of xval2                    " << endl;
        cout << "   -Yval=yval2      -> yval2 is the y value of the point              " << endl;
        cout << "   -Yerr=yerr2      -> yerr2 is the error of yval2                    " << endl;
        cout << "   -colP=index      -> color index of the measured point [default: index=1]" << endl;
        cout << "   -legP=legend     -> legend for the point [default: no legend]      " << endl;
        cout << "######################################################################" << endl;
        return 0;
    }

    //----  parameters which can be changed by command-line arguments  -----
    
    //double prob68 = 0.6827, prob95 = 0.9545, prob99 = 0.9973; //,0.9995,0.9999
    double prob68 = 0.68, prob95 = 0.95, prob99 = 0.99;
    //
    bool bOneDim = false, bCompat = false, bTwoDim = false;
    bool bPDF = false;
    bool bOrig = false, bOutputTxt = false, bContLines = false, bLeftLegend = false;
    bool bRescaleForMHl = false;
    int maxDig = 8, prec = 6;
    int nx = 100, ny = 20;
    double xval = -999.0, xerr = 0.0, x_low = 0.0, x_up = 0.0, y_low = 0.0, y_up = 0.0;
    TString addtext = "";
    double addtext_x = 0.0, addtext_y = 0.0;
    TString xlab = "", ylab = "";
    double yTitleOffset = 1.5;
    //
    double xval2 = -999.0, xerr2 = 0.0, yval2 = -999.0, yerr2 = 0.0;
    TString legP="";
    int colP = 1;
    //
    double legend_scale = 1.0;
    double leg_xmin = 0.63, leg_ymax = 0.88;
    TString legGauss = "";
    double prior_mean = 0.0, prior_sigma = 0.0;

    const int NumHist = 5;
    TObject* tobj[NumHist];
    TFile *datafile[NumHist];
    string plotname[NumHist], filename[NumHist];
    TString leg[NumHist];
    int smooth[NumHist], col68[NumHist], col95[NumHist];
    int lineStyle[NumHist], lineStyle68[NumHist], lineWidth[NumHist], lineColor[NumHist];
    double col68alpha[NumHist], col95alpha[NumHist];
    int fillStyle[NumHist], newNbins[NumHist], NumNewPoints[NumHist];
    bool bOnlyLine[NumHist], bOnly95[NumHist];
    for (int hist=0; hist<NumHist; hist++) {
        plotname[hist] = "";
        filename[hist] = "";
        leg[hist] = "";
        smooth[hist] = 0;
        lineStyle[hist] = 1;
        lineStyle68[hist] = 3;
        lineWidth[hist] = 1;
        lineColor[hist] = 1;
        col68[hist] = 1393;
        col95[hist] = 1392;
        col68alpha[hist] = 1.0;
        col95alpha[hist] = 1.0;
        fillStyle[hist] = 1001;
        newNbins[hist] = 100;
        NumNewPoints[hist] = 20;
        bOnlyLine[hist] = false;
        bOnly95[hist] = false;
    }

    //----------------------------------------------------------------------
    
    if (strncmp(argv[1], "--oneDim", 8) == 0) bOneDim = true;
    if (strncmp(argv[1], "--compat", 8) == 0) bCompat = true;
    if (strncmp(argv[1], "--twoDim", 8) == 0) bTwoDim = true;    

    if (strncmp(argv[2], "-rootfile1=", 11) == 0) {
        char str[100];
        sscanf(argv[2], "-rootfile1=%s", str);
        filename[0] = str;
    }        
    datafile[0] = new TFile(filename[0].c_str());
    if (argc == 3 && strncmp(argv[1], "info", 4) == 0) {
        datafile[0]->ls();
        return 0;
    }
    for (int hist=1; hist<NumHist; hist++)
        filename[hist] = filename[0];

    string outputFileName;   
    
    for (int i = 3; i < argc; i++) {
        char str[100];
        bool flag = true;
        
        if (strncmp(argv[i], "--orig", 6) == 0) bOrig = true;
        else if (strncmp(argv[i], "--pdf", 5) == 0) bPDF = true;
        else if (strncmp(argv[i], "--outputTxt", 11) == 0) bOutputTxt = true;
        else if (strncmp(argv[i], "--drawlines", 11) == 0) bContLines = true;
        else if (strncmp(argv[i], "--leftLegend", 12) == 0) bLeftLegend = true;
        else if (strncmp(argv[i], "--rescaleForMHl", 15) == 0) bRescaleForMHl = true;
        //
        else if (strncmp(argv[i], "-prob68=", 8) == 0)
            sscanf(argv[i], "-prob68=%lf", &prob68);
        else if (strncmp(argv[i], "-prob95=", 8) == 0)
            sscanf(argv[i], "-prob95=%lf", &prob95);
        else if (strncmp(argv[i], "-prob99=", 8) == 0)
            sscanf(argv[i], "-prob99=%lf", &prob99);
        //
        else if (strncmp(argv[i], "-output=", 8) == 0) {
            sscanf(argv[i], "-output=%s", str);
            outputFileName = str;
        }
        //
        else if (strncmp(argv[i], "-xlab=", 6) == 0) {
            sscanf(argv[i], "-xlab=%s", str);
            xlab = str;
        }
        //
        else if (strncmp(argv[i], "-ylab=", 6) == 0) {
            sscanf(argv[i], "-ylab=%s", str);
            ylab = str;
        }
        //
        else if (strncmp(argv[i], "-addtext=", 9) == 0) {
            sscanf(argv[i], "-addtext=%s", str);
            addtext = str;
        }
        else if (strncmp(argv[i], "-addtextAt=", 11) == 0)
            sscanf(argv[i], "-addtextAt=[%lf,%lf]", &addtext_x, &addtext_y);
        //
        else if (strncmp(argv[i], "-maxDigits=", 11) == 0)
            sscanf(argv[i], "-maxDigits=%d", &maxDig);
        //
        else if (strncmp(argv[i], "-precision=", 11) == 0)
            sscanf(argv[i], "-precision=%d", &prec);
        //
        else if (strncmp(argv[i], "-priorMean=", 11) == 0)
            sscanf(argv[i], "-priorMean=%lf", &prior_mean);
        else if (strncmp(argv[i], "-priorSigma=", 12) == 0)
            sscanf(argv[i], "-priorSigma=%lf", &prior_sigma);
        //
        else if (strncmp(argv[i], "-xrange=", 8) == 0) {
            TString stmp(argv[i] + 8);
            sscanf(stmp.Data(), "[%lf,%lf]", &x_low, &x_up);
        }
        //
        else if (strncmp(argv[i], "-range=", 7) == 0) {
            TString stmp(argv[i] + 7);
            sscanf(stmp.Data(), "[%lf,%lf]x[%lf,%lf]", &x_low, &x_up, &y_low, &y_up);
        }
        //
        else if (strncmp(argv[i], "-bins=", 6) == 0) {
            TString stmp(argv[i] + 6);
            sscanf(stmp.Data(), "[%d]x[%d]", &nx, &ny);
        }
        //
        else if (strncmp(argv[i], "-yTitleOffset=", 14) == 0)
            sscanf(argv[i], "-yTitleOffset=%lf", &yTitleOffset);
        //
        else if (strncmp(argv[i], "-val=", 5) == 0)
            sscanf(argv[i], "-val=%lf", &xval);
        else if (strncmp(argv[i], "-err=", 5) == 0)
            sscanf(argv[i], "-err=%lf", &xerr);
        //
        else if (strncmp(argv[i], "-Xval=", 6) == 0)
            sscanf(argv[i], "-Xval=%lf", &xval2);
        else if (strncmp(argv[i], "-Xerr=", 6) == 0)
            sscanf(argv[i], "-Xerr=%lf", &xerr2);
        else if (strncmp(argv[i], "-Yval=", 6) == 0)
            sscanf(argv[i], "-Yval=%lf", &yval2);
        else if (strncmp(argv[i], "-Yerr=", 6) == 0)
            sscanf(argv[i], "-Yerr=%lf", &yerr2);
        else if (strncmp(argv[i], "-colP=", 6) == 0)
            sscanf(argv[i], "-colP=%d", &colP);
        //
        else if (strncmp(argv[i], "-legScale=", 10) == 0)
            sscanf(argv[i], "-legScale=%lf", &legend_scale);
        else if (strncmp(argv[i], "-legXmin=", 9) == 0)
            sscanf(argv[i], "-legXmin=%lf", &leg_xmin);
        else if (strncmp(argv[i], "-legYmax=", 9) == 0)
            sscanf(argv[i], "-legYmax=%lf", &leg_ymax);
        else if (strncmp(argv[i], "-legP=", 6) == 0) {
            sscanf(argv[i], "-legP=%s", str);
            legP = str;
        } 
        else if (strncmp(argv[i], "-legGauss=", 10) == 0) {
            sscanf(argv[i], "-legGauss=%s", str);
            legGauss = str;
        }
        else
            flag = false;

        /* Options for multiple histograms */
        bool flag2 = true;
        for (int hist=1; hist<=NumHist; hist++) {
            stringstream ss;
            bool flagtmp = false;

            ss << "--onlyLine" << hist;
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                bOnlyLine[hist-1] = true;
                flagtmp = true;
            }
            ss.str("");

            ss << "--only95" << hist;
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                bOnly95[hist-1] = true;
                flagtmp = true;
            }
            ss.str("");
            
            ss << "-moreBins" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &newNbins[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-smooth" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &smooth[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-col68" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &col68[hist-1]);
                flagtmp = true;
            }
            ss.str("");
            
            ss << "-col95" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &col95[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-col68alpha" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%lf").c_str(), &col68alpha[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-col95alpha" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%lf").c_str(), &col95alpha[hist-1]);
                flagtmp = true;
            }
            ss.str("");
            
            ss << "-lineStyle" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &lineStyle[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-lineStyle68" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &lineStyle68[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-lineWidth" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &lineWidth[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-lineColor" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &lineColor[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-fill" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &fillStyle[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            ss << "-plot" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%s").c_str(), str);
                plotname[hist-1] = str;
                flagtmp = true;
            }
            ss.str("");

            ss << "-rootfile" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%s").c_str(), str);
                filename[hist-1] = str;
                flagtmp = true;
            }
            ss.str("");

            ss << "-leg" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%s").c_str(), str);
                leg[hist-1] = str;
                flagtmp = true;
            }
            ss.str("");

            ss << "-NumNewPoints" << hist << "=";
            if (strncmp(argv[i], ss.str().c_str(), ss.str().size()) == 0) {
                sscanf(argv[i], (ss.str() + "%d").c_str(), &NumNewPoints[hist-1]);
                flagtmp = true;
            }
            ss.str("");

            flag2 |= flagtmp;
        }

        if (!flag && !flag2) {
            cout << "Wrong option: " << argv[i] << endl;
            return 1;            
        }
    }

    // Errors
    if (!(bOneDim | bCompat | bTwoDim)) {
        cout << "Error: A mandatory option --oneDim, --compat or --twoDim is missing!" << endl;
        return 1;
    }
    if ( (bOneDim&bCompat) | (bOneDim&bTwoDim) | (bTwoDim&bCompat) ) {
        cout << "Error: Use only one of --oneDim, --compat and --twoDim!" << endl;
        return 1;
    }

    // output files
    if (outputFileName.compare("")==0)
        outputFileName = plotname[0];
    string epsFileName = outputFileName;
    if (!bPDF) epsFileName += ".eps";
    else epsFileName += ".pdf";
    string txtFileName = outputFileName + ".txt";    

    for (int n=0; n<NumHist; n++) {
        if ( plotname[n].compare("")!=0 ) {
            if (n!=0) datafile[n] = new TFile(filename[n].c_str());
            tobj[n] = datafile[n]->Get(plotname[n].c_str());
            if (tobj[n] == NULL) {
                cout << "Error: plot \"" << plotname[n] << "\" does not exist..." << endl;
                return 1;
            }
        }    
    }

    //----------------------------------------------------------------------
    // Colors

    BaseMacros::DefineNewColours();
    BaseMacros myMacros;

    // Fill colors
    TColor *FillCol68[NumHist], *FillCol95[NumHist];
    for (int i=0; i<NumHist; i++) {
        TColor *colTmp = gROOT->GetColor(col68[i]);
        if (gROOT->GetColor(col68[i]+10000) == NULL ) {
            FillCol68[i] = new TColor(col68[i]+10000, colTmp->GetRed(),
                                      colTmp->GetGreen(), colTmp->GetBlue());
            FillCol68[i]->SetAlpha(col68alpha[i]);
        } else
            FillCol68[i] = gROOT->GetColor(col68[i]+10000);
        TColor *colTmp2 = gROOT->GetColor(col95[i]);
        if (gROOT->GetColor(col95[i]+10000) == NULL ) {
            FillCol95[i] = new TColor(col95[i]+10000, colTmp2->GetRed(),
                                      colTmp2->GetGreen(), colTmp2->GetBlue());
            FillCol95[i]->SetAlpha(col95alpha[i]);
        } else
            FillCol95[i] = gROOT->GetColor(col95[i]+10000);
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

    //double CanvasWidth = 300.0;
    //double CanvasHeight = 300.0;
    //TCanvas TC("TC", "TC", CanvasWidth, CanvasHeight);
    //TC.SetWindowSize(CanvasWidth + (CanvasWidth - TC.GetWw()),
    //                 CanvasHeight + (CanvasHeight - TC.GetWh()));
    //gSystem->ProcessEvents();
    TCanvas TC("TC", "TC", 3);
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

    // New line styles
    gStyle->SetLineStyleString(11,"28 25");
    gStyle->SetLineStyleString(12,"35 25 13 25");
    
    TLegend *legend;
    int num_leg = 0;
    for (int n=0; n<NumHist; n++)
        if ( leg[n].CompareTo("")!= 0 ) num_leg++;
    if ( legP.CompareTo("")!= 0 ) num_leg++;
    if ( legGauss.CompareTo("")!= 0 ) num_leg++;
    double legend_ymin;
    if (num_leg == 1) legend_ymin = leg_ymax - 0.08;
    else if (num_leg == 2) legend_ymin = leg_ymax - 0.14;
    else if (num_leg == 3) legend_ymin = leg_ymax - 0.20;
    else if (num_leg == 4) legend_ymin = leg_ymax - 0.26;
    else if (num_leg == 5) legend_ymin = leg_ymax - 0.32;
    else legend_ymin = leg_ymax - 0.38;
    legend_ymin = (leg_ymax - 0.02)*(1.0 - legend_scale) + legend_ymin*legend_scale;
    if (bLeftLegend) leg_xmin -= 0.40; // default: 0.23
    legend = new TLegend(leg_xmin, legend_ymin, leg_xmin+0.12, leg_ymax);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.043*legend_scale);
    legend->SetMargin(0.7);

    //----------------------------------------------------------------------
    // 1-D histogram
    
    if (bOneDim) {
        SFH1D* SFHisto1D[NumHist];
        TH1D* plot_pt[NumHist];
        TH1D* hist[NumHist];

        hist[0] = (TH1D*) tobj[0]->Clone();   
        os << hist[0]->GetXaxis()->GetTitle() << " in " << plotname[0] << endl;

        SFHisto1D[0] = new SFH1D(*hist[0], prob68, prob95, prob99);
        SFHisto1D[0]->smoothHist(smooth[0]);
        SFHisto1D[0]->increaseNbins(newNbins[0]);
        SFHisto1D[0]->DrawAxes(xlab, ylab, maxDig, x_low, x_up, yTitleOffset); // draw the axes
        SFHisto1D[0]->Draw(lineStyle[0], lineWidth[0], lineColor[0], 
                           FillCol68[0]->GetNumber(), FillCol95[0]->GetNumber(),
                           fillStyle[0], bOnlyLine[0], bOrig);
        if (bOnlyLine[0]) 
            plot_pt[0] = SFHisto1D[0]->getNewHist();
        else 
            plot_pt[0] = SFHisto1D[0]->getNewHist68();
        
        // rescale (for mHl)
        if (bRescaleForMHl) {
            SFHisto1D[0]->getHistAxes()->Scale(10.0);
            SFHisto1D[0]->getNewHist()->Scale(10.0);
            SFHisto1D[0]->getNewHist68()->Scale(10.0);
            SFHisto1D[0]->getNewHist95()->Scale(10.0);
            //SFHisto1D[0]->getHistAxes()->GetXaxis()->SetRange(400,1400);
            leg[0] += " [x10]";
        }

        // superimpose other histograms
        for (int n=1; n<NumHist; n++) {
            if ( plotname[n].compare("")!=0 ) {
                hist[n] = (TH1D*) tobj[n]->Clone();
                SFHisto1D[n] = new SFH1D(*hist[n], prob68, prob95, prob99);
                SFHisto1D[n]->smoothHist(smooth[n]);
                SFHisto1D[n]->increaseNbins(newNbins[n]);
                SFHisto1D[n]->Draw(lineStyle[n], lineWidth[n], lineColor[n], 
                                   FillCol68[n]->GetNumber(), FillCol95[n]->GetNumber(),
                                   fillStyle[n], bOnlyLine[n], bOrig);
                if (bOnlyLine[n]) 
                    plot_pt[n] = SFHisto1D[n]->getNewHist();
                else 
                    plot_pt[n] = SFHisto1D[n]->getNewHist68();
                
                // rescale the y axis
                SFHisto1D[0]->RescaleYaxis(SFHisto1D[n]->getNewHist()->GetMaximum());
            } else 
                plot_pt[n] = NULL;
        }
            
        // superimpose a Gaussian (prior) function
        TF1* prior = NULL;
        if (prior_sigma != 0.0) {
            double xmin = SFHisto1D[0]->getHistAxes()->GetXaxis()->GetXmin();
            double xmax = SFHisto1D[0]->getHistAxes()->GetXaxis()->GetXmax();
            prior = new TF1("prior",
                    "1./sqrt(2.*TMath::Pi())/[1]* exp(- (x-[0])*(x-[0])/2./[1]/[1])",
                    xmin, xmax);
            prior->SetParameter(0, prior_mean);
            prior->SetParameter(1, prior_sigma);    
            prior->SetLineStyle(2);
            prior->SetLineWidth(4);
            prior->SetNpx(1000);
            prior->Draw("SAME");
            
            // rescale the y axis
            SFHisto1D[0]->RescaleYaxis(prior->GetMaximum());
        }

        // legends: Change the order if necessary. 
        string leg_opt;
        if (prior != NULL) legend->AddEntry(prior, myMacros.ConvertTitle(legGauss), "L");
        for (int n=1; n<NumHist; n++) {
            if (plot_pt[n] != NULL) {
                if (!bOnlyLine[n]) leg_opt = "F";
                else leg_opt = "L";
                legend->AddEntry(plot_pt[n], myMacros.ConvertTitle(leg[n]), leg_opt.c_str());
            }
        }
        if (!bOnlyLine[0]) leg_opt = "F";
        else leg_opt = "L";
        legend->AddEntry(plot_pt[0], myMacros.ConvertTitle(leg[0]), leg_opt.c_str());

        // output results to os
        SFHisto1D[0]->OutputResults(os, smooth[0], true);

        // output results for the original histogram before smoothing 
        // nor increasing the number of bins for comparison
        SFH1D orig1D(*hist[0], prob68, prob95, prob99);
        os << endl << "[Original histogram]" << endl;
        orig1D.OutputResults(os, 0, false);

        /* write the modified histograms to a root file */
        TFile *outRootFile;
        string outRootFileName = outputFileName + ".root";
        outRootFile = new TFile(outRootFileName.c_str(), "RECREATE");
        for (int n=0; n<NumHist; n++)
            if ( plotname[n].compare("")!=0 )
                SFHisto1D[n]->getNewHist()->Write(plotname[n].c_str());
        outRootFile->Close();
        cout << outRootFileName << " has been created." << endl;
        
    //----------------------------------------------------------------------
    // Compatibility plot
        
    } else if (bCompat) {
        TH1D* hist = (TH1D*) tobj[0]->Clone();
        os << hist->GetXaxis()->GetTitle() << " in " << plotname[0] << endl;
        os << "  Num of bins: " << nx << " x " << ny << endl;        
        
        Pull CompatPlot(*hist, nx, ny, x_low, x_up, y_low, y_up);
        CompatPlot.Draw(xlab, ylab, xval, xerr, maxDig, yTitleOffset);
        os << "  Pull: " << CompatPlot.calcPull(xval, xerr) << endl;
        
    //----------------------------------------------------------------------
    // 2-D histogram        

    } else if (bTwoDim) {    
        SFH2D* SFHisto2D[NumHist];
        TGraph* contour_pt[NumHist];
        TH2D* hist[NumHist];

        hist[0] = (TH2D*) tobj[0]->Clone();
        os << hist[0]->GetXaxis()->GetTitle() << " vs " << hist[0]->GetYaxis()->GetTitle() 
           << " in " << plotname[0] << endl;
        os << "  smooth: " << smooth[0] << " time(s)" << endl;
        
        SFHisto2D[0] = new SFH2D(*hist[0], os, prob68, prob95,
                                 x_low, x_up, y_low, y_up, NumNewPoints[0]);
        SFHisto2D[0]->smoothHist(smooth[0]);
        SFHisto2D[0]->Draw(xlab, ylab, lineWidth[0], lineColor[0],
                           FillCol68[0]->GetNumber(), FillCol95[0]->GetNumber(),
                           lineStyle[0], lineStyle68[0], fillStyle[0],
                           maxDig, bContLines, bOnly95[0], false, yTitleOffset);
        contour_pt[0] = SFHisto2D[0]->getContour();
        
        // superimpose other histograms
        for (int n=1; n<NumHist; n++) {
            if ( plotname[n].compare("")!=0 ) {
                os << "[Graph " << n+1 << "]" << endl;
                os << "  smooth: " << smooth[n] << " time(s)" << endl;
                hist[n] = (TH2D*) tobj[n]->Clone();
                SFHisto2D[n] = new SFH2D(*hist[n], os, prob68, prob95,
                                         x_low, x_up, y_low, y_up, NumNewPoints[n]);
                SFHisto2D[n]->smoothHist(smooth[n]);
                SFHisto2D[n]->Draw("", "", lineWidth[n], lineColor[n],
                                   FillCol68[n]->GetNumber(), FillCol95[n]->GetNumber(),
                                   lineStyle[n], lineStyle68[n], fillStyle[n], 
                                   maxDig, bContLines, bOnly95[n], true, yTitleOffset);
                contour_pt[n] = SFHisto2D[n]->getContour();
            } else 
                contour_pt[n] = NULL;
        }

        // draw a given point with error bars
        TGraphErrors* g1 = NULL;
        if (xval2 != -999.0 && yval2 != -999.0) {
            TLine *lx = new TLine();
            lx->SetLineWidth(3);
            lx->SetLineColor(colP);
            double zero = 0., err;
            err = xerr2;
            g1 = new TGraphErrors(1, &xval2, &yval2, &err, &zero);
            g1->SetLineWidth(3);
            g1->SetLineStyle(1);
            g1->SetLineColor(colP);
            g1->SetMarkerStyle(20);
            g1->SetMarkerSize(1);
            g1->SetMarkerColor(colP);
            g1->Draw("P");
            
            double xmin = hist[0]->GetXaxis()->GetXmin();
            double xmax = hist[0]->GetXaxis()->GetXmax();
            double ymin = hist[0]->GetYaxis()->GetXmin();
            double ymax = hist[0]->GetYaxis()->GetXmax();
            
            err = xerr2;
            double min_x = std::max(xval2 - err, xmin);
            double max_x = std::min(xval2 + err, xmax);
            lx->DrawLine(min_x, yval2, max_x, yval2);
            
            TLine *ly = new TLine();
            ly->SetLineWidth(3);
            ly->SetLineColor(colP);
            err = yerr2;
            TGraphErrors* g2 = new TGraphErrors(1, &xval2, &yval2, &zero, &err);
            g2->SetLineWidth(3);
            g2->SetLineStyle(1);
            g2->SetLineColor(colP);
            g2->SetMarkerStyle(20);
            g2->SetMarkerSize(1);
            g2->SetMarkerColor(colP);
            g2->Draw("P");
            
            double min_y = std::max(yval2 - err, ymin);
            double max_y = std::min(yval2 + err, ymax);
            ly->DrawLine(xval2, min_y, xval2, max_y);
        }
        
        // legends: Change the order if necessary. 
        string leg_opt;
        for (int n=NumHist-1; n>=0; n--) {
            if (contour_pt[n] != NULL && leg[n].CompareTo("")!= 0) {
                if (fillStyle[n]!=0) leg_opt = "F";
                else leg_opt = "L";
                legend->AddEntry(contour_pt[n], myMacros.ConvertTitle(leg[n]), leg_opt.c_str());
            }
        }
        if (g1 != NULL && legP.CompareTo("")!= 0) 
            legend->AddEntry(g1, myMacros.ConvertTitle(legP), "LP");
    } 

    //----------------------------------------------------------------------

    // additional text
    if ( addtext.CompareTo("")!= 0 ) {
        TLatex* TL = new TLatex();
        TL->SetTextSize(0.043);
        TL->SetTextFont(42);
        TL->SetTextAlign(11);
        TL->SetNDC(1);
        TL->DrawLatex(addtext_x, addtext_y, myMacros.ConvertTitle(addtext));
    }
       
    // draw the legend if not empty
    if (num_leg>0) legend->Draw("");
    
    gPad->RedrawAxis();
    TC.Print(epsFileName.c_str());

    /* close the output text file */
    if (bOutputTxt) {
        fout->close();          
        cout << txtFileName << " has been created." << endl;
    }

    return 0;
}

/** 
 * @}
 */


