/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TH1D.h>
#include <TStyle.h>
#include <TString.h>
#include "BaseMacros.h"

#ifndef SFH1D_H
#define	SFH1D_H

/**
 * @class SFH1D
 * @ingroup Macros
 * @brief A class for handling a 1-D histogram.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to draw the smallest intervals containing the 68% 
 * and 95% probability ranges. 
 */
class SFH1D : public BaseMacros {
public:

    /**
     * A constructor.
     * @param[in] hist A 1-D ROOT histogram 
     * @param[in] prob68_in The probability of the 68% interval. 
     * @param[in] prob95_in The probability of the 95% interval. 
     */
    SFH1D(TH1D& hist, const double prob68_in=0.68, const double prob95_in=0.95);
 
    /**
     * @return Empty histogram for the axes.
     */
    TH1D* getHistAxes() const 
    {
        return HistAxes;
    }

    /**
     * @return The modified 1-D histogram. 
     */
    TH1D* getNewHist() const 
    {
        return newHist;
    }
    
    /**
     * @return A histogram containing only the 68% interval. 
     */
    TH1D* getNewHist68() const 
    {
        return newHist68;
    }
    
    /**
     * @return A histogram containing only the 95% interval. 
     */
    TH1D* getNewHist95() const 
    {
        return newHist95;
    }
    
    /**
     * @return The minimum of the 68% interval.
     */ 
    double getXmin68() const 
    {
        return xmin68;
    }

    /**
     * @return The maximum of the 68% interval.
     */
    double getXmax68() const 
    {
        return xmax68;
    }

    /**
     * @return The minimum of the 95% interval.
     */
    double getXmin95() const 
    {
        return xmin95;
    }

    /**
     * @return The maximum of the 95% interval.
     */
    double getXmax95() const 
    {
        return xmax95;
    }
    
    /**
     * @return The local mode of the modified histogram. 
     */
    double getLocalMode() const
    {
        return localMode;        
    }    
    
    /**
     * Smooth the bin contents of th histogram by using TH1::Smooth(). 
     * @param smooth
     */
    void smoothHist(const int smooth=0);
    
    /**
     * Increase the number of bins, where the contents of new bins are fitted 
     * with quadratic functions. 
     * @param[in] newNbin The number of new bins [Default: 1000]. 
     */
    void increaseNbins(const int newNbin=1000);

    /**
     * Draw the axes. 
     * @param[in] xlab The label of the x axis. 
     * @param[in] ylab The label of the y axis. 
     * @param[in] maxDigits The maximum digits of axis labels. 
     * @param[in] x_low
     * @param[in] x_up
     * @param[in] YTitleOffset 
     */
    void DrawAxes(const TString xlab, const TString ylab, 
                  const int maxDigits, 
                  const double x_low = 0.0, const double x_up = 0.0, 
                  const double YTitleOffset=1.5);
    
    /**
     * Draw the modified histogram. 
     * @param[in] lineStyle The index of the line style. 
     * @param[in] lineWidth The line width. 
     * @param[in] lineCorlo The index of the line color. 
     * @param[in] col68 The color index for the 68% interval. 
     * @param[in] col95 The color index for the 95% interval. 
     * @param[in] fillStyle The index of the fill area style. 
     * @param[in] bOnlyLine 
     * @param[in] bOrigHist A flag controlling if the original histogram is superimposed. 
     */
    void Draw(const int lineStyle, const int lineWidth, const int lineColor, 
              const int col68, const int col95, const int fillStyle,
              const bool bOnlyLine = false, const bool bOrigHist = false);

    /**
     * Rescale the y axis.
     * @param[in] max_another The maximum of another histogram to be superimposed. 
     */    
    void RescaleYaxis(double max_another);
    
    /**
     * Output results. 
     * @param[in] os An output stream. 
     * @param[in] smooth 
     * @param[in] WasDrawed
     */    
    void OutputResults(ostream& os, const int smooth, const bool WasDrawed = false) const;
    
    
private:

    /**
     * Probability of the 68% interval.
     */
    const double prob68;
    
    /**
     * Probability of the 95% interval. 
     */
    const double prob95;

    /**
     * The original 1-D histogram.
     */
    TH1D& origHist;

    /**
     * The name of the original histogram.
     */
    const std::string origName;
    
    /*
     * Empty histogram for the axes.
     */
    TH1D* HistAxes;

    /**
     * A modified histogram.  
     */
    TH1D* newHist;
    
    /**
     * A histogram containing only the 68% interval of newHist. 
     */
    TH1D *newHist68;
    
    /**
     * A histogram containing only the 95% interval of newHist. 
     */
    TH1D *newHist95;    

    /**
     * The minimum of the 68% interval.
     */
    double xmin68;
    
    /**
     * The maximum of the 68% interval.
     */
    double xmax68;
    
    /**
     * The minimum of the 95% interval.
     */
    double xmin95;
    
    /**
     * The maximum of the 95% interval.
     */
    double xmax95;
    
    /**
     * The local mode of the modified histogram. 
     */
    double localMode;
    
    /**
     * Compute @f$y(x)=ax^2+bx+c@f$ for a given x, where a, b and c are fixed 
     * from the pairs of xx and yy.
     * @param[in] x A point to be interested. 
     * @param[in] xx An array of x with 3 elements. 
     * @param[in] yy An array of y(x) with 3 elements. 
     * @return The value at the given point x. 
     */
    double quadra(const double x, const double xx[3], const double yy[3]) const;   

    /**
     * Create a new histogram containing the smallest interval at a given level. 
     * @param[out] min The minimum of the interval. 
     * @param[out] max The maximum of the interval. 
     * @param[in] level A probability level. 
     * @return A 1-D histogram containing the smallest interval at a given level. 
     */
    TH1D* HistInterval(double &min, double &max, const double level) const;

    /**
     * Compute the 68% and 95% intervals with  
     */
    void computeIntervals();
    
};

#endif	/* SFH1D_H */

