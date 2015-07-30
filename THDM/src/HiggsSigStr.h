/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSSIGSTR_H
#define	HIGGSSIGSTR_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class HiggsSigStr
 * @ingroup THDM 
 * @brief An observable class for the 2HDM Higgs signal strengths.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute ten 2HDM Higgs signal strengths.
 */
class HiggsSigStr : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   HiggsSigStr(const StandardModel& SM_i, int obsFlag);
     
    /**
     * @brief Ten 2HDM Higgs signal strengths.
     * @return
     */
    double computeThValue();

    
    private:
        const THDM * myTHDM;
        int obs;

    protected:
        gslpp::complex f_func(const double x) const;
        gslpp::complex g_func(const double x) const;
        gslpp::complex Int1(const double tau, const double lambda) const;
        gslpp::complex Int2(const double tau, const double lambda) const;
};

#endif	/* HIGGSSIGSTR_H */
