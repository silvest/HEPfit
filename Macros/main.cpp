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
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -prob68=prob     -> probability for the first interval [default: prob=0.68]" << endl;
        cout << "   -prob95=prob     -> probability for the second interval [default: prob=0.95]"<< endl;
        cout << "   -xlab=namex      -> x label                                        " << endl;
        cout << "   -ylab=namey      -> y label                                        " << endl;
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -xrange=\"[xmin,xmax]\" -> define the graph range                  " << endl;
        cout << "   -legScale=scale  -> scale factor for the legend [default: scale=1] " << endl;
        cout << "   -legXmin=xmin    -> xmin for the position of the legend [default: xmin=0.63]" << endl;
        cout << "   -legYmax=ymax    -> ymax for the position of the legend [default: ymax=0.88]" << endl;
        cout << "   *** options for histograms (N=1,2,3,4) ***                         " << endl;
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
        cout << "   -fillN=index     -> index of the fill area style [default: index=1001]" << endl;
        cout << "   -legN=legend     -> legend for plotN [default: no legend]          " << endl;
        cout << "   [e.g. -plot2=Mw -roogfile2=test.root -smooth2=3 -col682=2 -col952=5]" << endl;
        cout << "   *** superimpose a Gaussian function ***                            " << endl;
        cout << "   -priorMean=mean  -> mean value for the Gaussian function           " << endl;
        cout << "   -priorSigma=sig  -> standard deviation for the Gaussian function   " << endl;
        cout << "   -legGauss=legend -> legend for the Gaussian function [default: no legend]" << endl;
        cout << "                                                                      " << endl;
        cout << " Optional parameters for compatibility plots:                         " << endl;
        cout << "   --compat         -> Compatibility plot (mandatory)                 " << endl;
        cout << "   --outputTxt      -> output results to a text file                  " << endl;
        cout << "   -output=filename -> the base name of output eps and text files (without extension)" << endl;
        cout << "                       [default: plotname]                            " << endl;
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -prob68=prob     -> probability for the first interval [default: prob=0.68]" << endl;
        cout << "   -prob95=prob     -> probability for the second interval [default: prob=0.95]"<< endl;
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range       " << endl;
        cout << "   -bins=\"[xbins]x[ybins]\"  -> define the graph binning [default: 100x20]" << endl;
        cout << "   -xlab=namex      -> x label                                        " << endl;
        cout << "   -ylab=namey      -> y label                                        " << endl;
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
        cout << "   -maxDigits=n     -> set max digits in the axis labels [default: n=8]" << endl;
        cout << "   -precision=n     -> precision of values in the std output [default: n=6]" << endl;
        cout << "   -addtext=text    -> attach additional information                  " << endl;
        cout << "   -addtextAt=\"[x,y]\" -> position of the text                       " << endl;
        cout << "   -range=\"[xmin,xmax]x[ymin,ymax]\" -> define the graph range       " << endl;
        cout << "   -xlab=namex      -> x label                                        " << endl;
        cout << "   -ylab=namey      -> y label                                        " << endl;
        cout << "   -legScale=scale  -> scale factor for the legend [default: scale=1.0]" << endl;
        cout << "   -legXmin=xmin    -> xmin for the position of the legend [default: xmin=0.63]" << endl;
        cout << "   -legYmax=ymax    -> ymax for the position of the legend [default: ymax=0.88]" << endl;
        cout << "   *** options for histogram (N=1,2,3,4) ***                          " << endl;
        cout << "   --only95N        -> draw only the 95% contour                      " << endl;
        cout << "   -plotN=name      -> name of the histogram                          " << endl;
        cout << "   -rootfileN=filename -> rootfile name for plotN (with extension)    " << endl;
        cout << "         [For N=2,3, the same as the rootfile for plot1 by default]   " << endl;
        cout << "   -smoothN=ntime   -> iterative smoothing for plotN [default: ntime=0]" << endl;
        cout << "   -col68N=index    -> color index of the 68% interval for plotN [default: index=1393]" << endl;
        cout << "   -col95N=index    -> color index of the 95% interval for plotN [default: index=1392]" << endl;
        cout << "   -lineN=index     -> index of the line style [default: index=1]     " << endl;
        cout << "   -fillN=index     -> index of the fill area style [default: index=1001]" << endl;
        cout << "   -legN=legend     -> legend for plotN [default: no legend]          " << endl;
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
    
    //double prob68 = 0.6827, prob95 = 0.9545;
    double prob68 = 0.68, prob95 = 0.95;
    //
    bool bOneDim = false, bCompat = false, bTwoDim = false;
    bool bOrig = false, bOutputTxt = false, bContLines = false, bLeftLegend = false;
    int maxDig = 8, prec = 6;
    int nx = 100, ny = 20;
    double xval = -999.0, xerr = 0.0, x_low = 0.0, x_up = 0.0, y_low = 0.0, y_up = 0.0;
    TString addtext = "";
    double addtext_x = 0.0, addtext_y = 0.0;
    TString xlab = "", ylab = "";
    //
    double xval2 = -999.0, xerr2 = 0.0, yval2 = -999.0, yerr2 = 0.0;
    TString legP="";
    int colP = 1;
    //
    TString legGauss = "";
    double prior_mean = 0.0, prior_sigma = 0.0;
    //
    const int NumHist = 4;   
    TObject* tobj[NumHist];
    TFile *datafile[NumHist];
    string plotname[NumHist] = {"", "", "", ""};
    string filename[NumHist] = {"", "", "", ""};
    TString leg[NumHist] = {"", "", "", ""};
    int smooth[NumHist] = {0, 0, 0, 0};
    int lineStyle[NumHist] = {1, 1, 1, 1};
    int lineWidth[NumHist] = {1, 1, 1, 1};
    int lineColor[NumHist] = {1, 1, 1, 1};
    int col68[NumHist] = {1393, 1393, 1393, 1393};
    int col95[NumHist] = {1392, 1392, 1392, 1392};
    int fillStyle[NumHist] = {1001, 1001, 1001, 1001};
    int newNbins[NumHist] = {100, 100, 100, 100};
    bool bOnlyLine[NumHist] = {false, false, false, false};
    bool bOnly95[NumHist] = {false, false, false, false};
    
    double legend_scale = 1.0;
    double leg_xmin = 0.63, leg_ymax = 0.88;

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
    filename[1] = filename[0];
    filename[2] = filename[0];
    filename[3] = filename[0];

    string outputFileName;   
    
    for (int i = 3; i < argc; i++) {
        char str[100];
        if (strncmp(argv[i], "--orig", 6) == 0) bOrig = true;
        else if (strncmp(argv[i], "--outputTxt", 11) == 0) bOutputTxt = true;        
        else if (strncmp(argv[i], "--drawlines", 11) == 0) bContLines = true;        
        else if (strncmp(argv[i], "--leftLegend", 12) == 0) bLeftLegend = true;        
        else if (strncmp(argv[i], "--onlyLine1", 11) == 0) bOnlyLine[0] = true;
        else if (strncmp(argv[i], "--onlyLine2", 11) == 0) bOnlyLine[1] = true;
        else if (strncmp(argv[i], "--onlyLine3", 11) == 0) bOnlyLine[2] = true;
        else if (strncmp(argv[i], "--onlyLine4", 11) == 0) bOnlyLine[3] = true;
        else if (strncmp(argv[i], "--only951", 9) == 0) bOnly95[0] = true;
        else if (strncmp(argv[i], "--only952", 9) == 0) bOnly95[1] = true; 
        else if (strncmp(argv[i], "--only953", 9) == 0) bOnly95[2] = true;
        else if (strncmp(argv[i], "--only954", 9) == 0) bOnly95[3] = true;

        else if (strncmp(argv[i], "-prob68=", 8) == 0) 
            sscanf(argv[i], "-prob68=%lf", &prob68);
        else if (strncmp(argv[i], "-prob95=", 8) == 0) 
            sscanf(argv[i], "-prob95=%lf", &prob95);
        
        else if (strncmp(argv[i], "-output=", 8) == 0) {
            sscanf(argv[i], "-output=%s", str);
            outputFileName = str;
        }

        else if (strncmp(argv[i], "-xlab=", 6) == 0) {
            sscanf(argv[i], "-xlab=%s", str);
            xlab = str;
        }
        
        else if (strncmp(argv[i], "-ylab=", 6) == 0) {
            sscanf(argv[i], "-ylab=%s", str);
            ylab = str;
        }        
        
        else if (strncmp(argv[i], "-addtext=", 9) == 0) {
            sscanf(argv[i], "-addtext=%s", str);
            addtext = str;
        }
        else if (strncmp(argv[i], "-addtextAt=", 11) == 0)
            sscanf(argv[i], "-addtextAt=[%lf,%lf]", &addtext_x, &addtext_y);
        
        else if (strncmp(argv[i], "-maxDigits=", 11) == 0) 
            sscanf(argv[i], "-maxDigits=%d", &maxDig);

        else if (strncmp(argv[i], "-precision=", 11) == 0) 
            sscanf(argv[i], "-precision=%d", &prec);
        
        else if (strncmp(argv[i], "-priorMean=", 11) == 0) 
            sscanf(argv[i], "-priorMean=%lf", &prior_mean);
        else if (strncmp(argv[i], "-priorSigma=", 12) == 0) 
            sscanf(argv[i], "-priorSigma=%lf", &prior_sigma);

        else if (strncmp(argv[i], "-xrange=", 8) == 0) {
            TString stmp(argv[i] + 8);
            sscanf(stmp.Data(), "[%lf,%lf]", &x_low, &x_up);
        }

        else if (strncmp(argv[i], "-range=", 7) == 0) {
            TString stmp(argv[i] + 7);
            sscanf(stmp.Data(), "[%lf,%lf]x[%lf,%lf]", &x_low, &x_up, &y_low, &y_up);
        }

        else if (strncmp(argv[i], "-bins=", 6) == 0) {
            TString stmp(argv[i] + 6);
            sscanf(stmp.Data(), "[%d]x[%d]", &nx, &ny);
        }

        else if (strncmp(argv[i], "-val=", 5) == 0) 
            sscanf(argv[i], "-val=%lf", &xval);
        else if (strncmp(argv[i], "-err=", 5) == 0) 
            sscanf(argv[i], "-err=%lf", &xerr);
        
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

        else if (strncmp(argv[i], "-moreBins1=", 11) == 0)
            sscanf(argv[i], "-moreBins1=%d", &newNbins[0]);
        else if (strncmp(argv[i], "-moreBins2=", 11) == 0)
            sscanf(argv[i], "-moreBins2=%d", &newNbins[1]);
        else if (strncmp(argv[i], "-moreBins3=", 11) == 0)
            sscanf(argv[i], "-moreBins3=%d", &newNbins[2]);
        else if (strncmp(argv[i], "-moreBins4=", 11) == 0)
            sscanf(argv[i], "-moreBins4=%d", &newNbins[3]);
        
        else if (strncmp(argv[i], "-smooth1=", 9) == 0)
            sscanf(argv[i], "-smooth1=%d", &smooth[0]);
        else if (strncmp(argv[i], "-smooth2=", 9) == 0) 
            sscanf(argv[i], "-smooth2=%d", &smooth[1]);
        else if (strncmp(argv[i], "-smooth3=", 9) == 0) 
            sscanf(argv[i], "-smooth3=%d", &smooth[2]);
        else if (strncmp(argv[i], "-smooth4=", 9) == 0) 
            sscanf(argv[i], "-smooth4=%d", &smooth[3]);

        else if (strncmp(argv[i], "-col681=", 8) == 0) 
            sscanf(argv[i], "-col681=%d", &col68[0]);
        else if (strncmp(argv[i], "-col951=", 8) == 0) 
            sscanf(argv[i], "-col951=%d", &col95[0]);
        else if (strncmp(argv[i], "-col682=", 8) == 0) 
            sscanf(argv[i], "-col682=%d", &col68[1]);
        else if (strncmp(argv[i], "-col952=", 8) == 0) 
            sscanf(argv[i], "-col952=%d", &col95[1]);
        else if (strncmp(argv[i], "-col683=", 8) == 0) 
            sscanf(argv[i], "-col683=%d", &col68[2]);
        else if (strncmp(argv[i], "-col953=", 8) == 0) 
            sscanf(argv[i], "-col953=%d", &col95[2]);
        else if (strncmp(argv[i], "-col684=", 8) == 0) 
            sscanf(argv[i], "-col684=%d", &col68[3]);
        else if (strncmp(argv[i], "-col954=", 8) == 0) 
            sscanf(argv[i], "-col954=%d", &col95[3]);

        else if (strncmp(argv[i], "-lineStyle1=", 12) == 0) 
            sscanf(argv[i], "-lineStyle1=%d", &lineStyle[0]);
        else if (strncmp(argv[i], "-lineStyle2=", 12) == 0) 
            sscanf(argv[i], "-lineStyle2=%d", &lineStyle[1]);
        else if (strncmp(argv[i], "-lineStyle3=", 12) == 0) 
            sscanf(argv[i], "-lineStyle3=%d", &lineStyle[2]);
        else if (strncmp(argv[i], "-lineStyle4=", 12) == 0) 
            sscanf(argv[i], "-lineStyle4=%d", &lineStyle[3]);
        
        else if (strncmp(argv[i], "-lineWidth1=", 12) == 0) 
            sscanf(argv[i], "-lineWidth1=%d", &lineWidth[0]);
        else if (strncmp(argv[i], "-lineWidth2=", 12) == 0) 
            sscanf(argv[i], "-lineWidth2=%d", &lineWidth[1]);
        else if (strncmp(argv[i], "-lineWidth3=", 12) == 0) 
            sscanf(argv[i], "-lineWidth3=%d", &lineWidth[2]);
        else if (strncmp(argv[i], "-lineWidth4=", 12) == 0) 
            sscanf(argv[i], "-lineWidth4=%d", &lineWidth[3]);

        else if (strncmp(argv[i], "-lineColor1=", 12) == 0) 
            sscanf(argv[i], "-lineColor1=%d", &lineColor[0]);
        else if (strncmp(argv[i], "-lineColor2=", 12) == 0) 
            sscanf(argv[i], "-lineColor2=%d", &lineColor[1]);
        else if (strncmp(argv[i], "-lineColor3=", 12) == 0) 
            sscanf(argv[i], "-lineColor3=%d", &lineColor[2]);
        else if (strncmp(argv[i], "-lineColor4=", 12) == 0) 
            sscanf(argv[i], "-lineColor4=%d", &lineColor[3]);

        else if (strncmp(argv[i], "-fill1=", 7) == 0) 
            sscanf(argv[i], "-fill1=%d", &fillStyle[0]);
        else if (strncmp(argv[i], "-fill2=", 7) == 0) 
            sscanf(argv[i], "-fill2=%d", &fillStyle[1]);
        else if (strncmp(argv[i], "-fill3=", 7) == 0) 
            sscanf(argv[i], "-fill3=%d", &fillStyle[2]);
        else if (strncmp(argv[i], "-fill4=", 7) == 0) 
            sscanf(argv[i], "-fill4=%d", &fillStyle[3]);
        
        else if (strncmp(argv[i], "-plot1=", 7) == 0) {
            sscanf(argv[i], "-plot1=%s", str);
            plotname[0] = str;
        }
        else if (strncmp(argv[i], "-plot2=", 7) == 0) {
            sscanf(argv[i], "-plot2=%s", str);
            plotname[1] = str;
        }
        else if (strncmp(argv[i], "-plot3=", 7) == 0) {
            sscanf(argv[i], "-plot3=%s", str);
            plotname[2] = str;
        }
        else if (strncmp(argv[i], "-plot4=", 7) == 0) {
            sscanf(argv[i], "-plot4=%s", str);
            plotname[3] = str;
        }
        
        else if (strncmp(argv[i], "-rootfile2=", 11) == 0) {
            sscanf(argv[i], "-rootfile2=%s", str);
            filename[1] = str;
        }        
        else if (strncmp(argv[i], "-rootfile3=", 11) == 0) {
            sscanf(argv[i], "-rootfile3=%s", str);
            filename[2] = str;
        }        
        else if (strncmp(argv[i], "-rootfile4=", 11) == 0) {
            sscanf(argv[i], "-rootfile4=%s", str);
            filename[3] = str;
        }        

        else if (strncmp(argv[i], "-legScale=", 10) == 0) 
            sscanf(argv[i], "-legScale=%lf", &legend_scale);
        else if (strncmp(argv[i], "-legXmin=", 9) == 0) 
            sscanf(argv[i], "-legXmin=%lf", &leg_xmin);
        else if (strncmp(argv[i], "-legYmax=", 9) == 0) 
            sscanf(argv[i], "-legYmax=%lf", &leg_ymax);
        else if (strncmp(argv[i], "-leg1=", 6) == 0) {
            sscanf(argv[i], "-leg1=%s", str);
            leg[0] = str;
        }
        else if (strncmp(argv[i], "-leg2=", 6) == 0) {
            sscanf(argv[i], "-leg2=%s", str);
            leg[1] = str;
        }
        else if (strncmp(argv[i], "-leg3=", 6) == 0) {
            sscanf(argv[i], "-leg3=%s", str);
            leg[2] = str;
        }
        else if (strncmp(argv[i], "-leg4=", 6) == 0) {
            sscanf(argv[i], "-leg4=%s", str);
            leg[3] = str;
        }
        else if (strncmp(argv[i], "-legGauss=", 10) == 0) {
            sscanf(argv[i], "-legGauss=%s", str);
            legGauss = str;
        }
        else if (strncmp(argv[i], "-legP=", 6) == 0) {
            sscanf(argv[i], "-legP=%s", str);
            legP = str;
        }
        
        else {
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
    string epsFileName = outputFileName + ".eps";
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
    BaseMacros myMacros;
    
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

        SFHisto1D[0] = new SFH1D(*hist[0], prob68, prob95);
        SFHisto1D[0]->smoothHist(smooth[0]);
        SFHisto1D[0]->increaseNbins(newNbins[0]);
        SFHisto1D[0]->DrawAxes(xlab, ylab, maxDig, x_low, x_up); // draw the axes
        SFHisto1D[0]->Draw(lineStyle[0], lineWidth[0], lineColor[0], 
                           col68[0], col95[0], fillStyle[0], bOnlyLine[0], bOrig); 
        if (bOnlyLine[0]) 
            plot_pt[0] = SFHisto1D[0]->getNewHist();
        else 
            plot_pt[0] = SFHisto1D[0]->getNewHist68();
        
        // rescale (for mHl)
        //SFHisto1D[0]->getHistAxes()->Scale(10.0);
        //SFHisto1D[0]->getNewHist()->Scale(10.0);
        //SFHisto1D[0]->getNewHist68()->Scale(10.0);
        //SFHisto1D[0]->getNewHist95()->Scale(10.0);
        //SFHisto1D[0]->getHistAxes()->GetXaxis()->SetRange(400,1400);
        //leg[0] += " [x10]";
            
        // superimpose other histograms
        for (int n=1; n<NumHist; n++) {
            if ( plotname[n].compare("")!=0 ) {
                hist[n] = (TH1D*) tobj[n]->Clone();
                SFHisto1D[n] = new SFH1D(*hist[n], prob68, prob95);
                SFHisto1D[n]->smoothHist(smooth[n]);
                SFHisto1D[n]->increaseNbins(newNbins[n]);
                SFHisto1D[n]->Draw(lineStyle[n], lineWidth[n], lineColor[n], 
                                   col68[n], col95[n], fillStyle[n], bOnlyLine[n], bOrig);    
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
        if (plot_pt[1] != NULL) {
            if (!bOnlyLine[1]) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(plot_pt[1], myMacros.ConvertTitle(leg[1]), leg_opt.c_str());
        }
        if (plot_pt[2] != NULL) {
            if (!bOnlyLine[2]) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(plot_pt[2], myMacros.ConvertTitle(leg[2]), leg_opt.c_str());
        }
        if (plot_pt[3] != NULL) {
            if (!bOnlyLine[3]) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(plot_pt[3], myMacros.ConvertTitle(leg[3]), leg_opt.c_str());
        }
        if (!bOnlyLine[0]) leg_opt = "F";
        else leg_opt = "L";
        legend->AddEntry(plot_pt[0], myMacros.ConvertTitle(leg[0]), leg_opt.c_str());
        
        // output results to os
        SFHisto1D[0]->OutputResults(os, smooth[0], true);

        // output results for the original histogram before smoothing 
        // nor increasing the number of bins for comparison
        SFH1D orig1D(*hist[0], prob68, prob95);
        os << endl << "[Original histogram]" << endl;
        orig1D.OutputResults(os, 0, false);
 
    //----------------------------------------------------------------------
    // Compatibility plot
        
    } else if (bCompat) {
        TH1D* hist = (TH1D*) tobj[0]->Clone();
        os << hist->GetXaxis()->GetTitle() << " in " << plotname[0] << endl;
        os << "  Num of bins: " << nx << " x " << ny << endl;
        
        Pull CompatPlot(*hist, nx, ny, x_low, x_up, y_low, y_up);
        CompatPlot.Draw(xlab, ylab, xval, xerr, maxDig);
        
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
        
        SFHisto2D[0] = new SFH2D(*hist[0], os, prob68, prob95, x_low, x_up, y_low, y_up);
        SFHisto2D[0]->smoothHist(smooth[0]);
        SFHisto2D[0]->draw(xlab, ylab, col68[0], col95[0], lineStyle[0], fillStyle[0], 
                           maxDig, bContLines, bOnly95[0], false);
        contour_pt[0] = SFHisto2D[0]->getContour();
        
        // superimpose other histograms
        for (int n=1; n<NumHist; n++) {
            if ( plotname[n].compare("")!=0 ) {
                os << "[Graph " << n+1 << "]" << endl;
                os << "  smooth: " << smooth[n] << " time(s)" << endl;
                hist[n] = (TH2D*) tobj[n]->Clone();
                SFHisto2D[n] = new SFH2D(*hist[n], os, prob68, prob95);
                SFHisto2D[n]->smoothHist(smooth[n]);
                SFHisto2D[n]->draw("", "", col68[n], col95[n], 
                                   lineStyle[n], fillStyle[n], 
                                   maxDig, bContLines, bOnly95[n], true);
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
            double zero = 0, err;
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
        if (contour_pt[3] != NULL && leg[3].CompareTo("")!= 0) {
            if (fillStyle[3]!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour_pt[3], myMacros.ConvertTitle(leg[3]), leg_opt.c_str());
        }
        if (contour_pt[2] != NULL && leg[2].CompareTo("")!= 0) {
            if (fillStyle[2]!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour_pt[2], myMacros.ConvertTitle(leg[2]), leg_opt.c_str()); 
        }
        if (contour_pt[1] != NULL && leg[1].CompareTo("")!= 0) {
            if (fillStyle[1]!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour_pt[1], myMacros.ConvertTitle(leg[1]), leg_opt.c_str());
        }
        if (leg[0].CompareTo("")!= 0) {
            if (fillStyle[0]!=0) leg_opt = "F";
            else leg_opt = "L";
            legend->AddEntry(contour_pt[0], myMacros.ConvertTitle(leg[0]), leg_opt.c_str());
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
    
    if (bOutputTxt) {
        fout->close();          
        cout << txtFileName << " has been created." << endl;
    }

    return 0;
}

/** 
 * @}
 */


