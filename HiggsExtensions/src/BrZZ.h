/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BrZZ_H
#define	BrZZ_H
#include <ThObservable.h>

/**
 * @class BrZZ
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to ZZ@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to ZZ@f]
 * in the current model and in the Standard Model
 */

class BrZZ : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrZZ(const HiggsExtensionModel& HESM_i) {
        if (HESM_i.ModelName().compare(0, 5, "Higgs") == 0)
            HESM = static_cast<HiggsExtensionModel&> (HESM_i);
        else
            throw std::runtime_error("ERROR: the BrZZ constructor can only be used with a HiggsExtensionModel reference, "
                + "while I got " + HESM_i.ModelName() + " as argument");

    };
    BrZZ(const BrZZ& orig);
    virtual ~BrZZ();

    /**
     * method to compute the the ratio of the @f[BR(H\to ZZ@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKZ() * HESM.computeKZ() / HESM.computeGTotalRatio();
    }
private:
    HiggsExtensionModel& HESM;
};

#endif	/* BrZZ_H */
