/*
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMLOWMASS_H
#define GENERALTHDMLOWMASS_H

#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class GeneralTHDMLowMass
 * @ingroup GeneralTHDM 
 * @brief Base class for low mass Higgs direct search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theory over experiment ratios of
 * several phenomena that include an extra scalar lighter than the SM Higgs
 */

/**
 * @class Hobs_pp_h_phi3phi3_mumutautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to \mu\mu\tau\tau@f$.
 */
class Hobs_pp_h_phi3phi3_mumutautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_mumutautau_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_mumutautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM/SM}}_{pp\to h_{125}}\cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \mu\mu) \cdot BR^{\text{GTHDM}}(A\to \tau\tau)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

#endif /* GENERALTHDMLOWMASS_H */