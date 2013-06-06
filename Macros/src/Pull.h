/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TLine.h>
#include <TMarker.h>
#include "BaseMacros.h"

#ifndef PULL_H
#define	PULL_H

/**
 * @class Pull
 * @ingroup Macros
 * @brief A class for drawing a compatibility plot. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Pull : public BaseMacros {
public:
    
    /**
     * A constructor. 
     * @param[in] hist A 1-D ROOT histogram 
     * @param[in] nbinX
     * @param[in] nbinY
     * @param[in] xLow
     * @param[in] xUp
     * @param[in] yLow
     * @param[in] yUp
     */
    Pull(TH1D& hist, const int nbinX, const int nbinY, 
         const double xLow, const double xUp, 
         const double yLow, const double yUp); 

    /**
     * Draw a compatibility plot. 
     * @param[in] xlab The label of the x axis. 
     * @param[in] ylab The label of the y axis. 
     * @param[in] xval The measured point to be superimposed. 
     * @param[in] xerr The error of xval. 
     * @param[in] maxDigits The maximum digits of axis labels. 
     * @param[in] YTitleOffset 
     */
    void Draw(const TString xlab, const TString ylab, 
              const double xval, const double xerr, const int maxDigits, 
              const double YTitleOffset=1.5);

    double calcPull(const double mean, const double sigma, const bool lowStat=false);



private:

    int nx;
    int ny;
    double x_low;
    double x_up;
    double y_low;
    double y_up;

    /**
     * A clone of a given 1-D histogram.   
     */
    TH1D* copyHist;

    TH2D* CompatPlot;
    
    TLine* lx;
    
    TLatex* tText[6];

    TMarker* ExpData;
    
    double meanTmp;
    
    double sigmaTmp;
    
    double integrand(const double* x) const;

    /**
     * @param[in] x
     * @return The content of the bin for a given x, where those of the 
     * underflow and overflow bins are set to 0. 
     */
    double fhisto(const double x) const;    
    
    void makeCompatPlot();
    
};

#endif	/* PULL_H */

