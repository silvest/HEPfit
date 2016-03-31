/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef OUTPUTSLHAFROMFH_H
#define	OUTPUTSLHAFROMFH_H

#include <ThObservable.h>


#include "FeynHiggsWrapper.h"

/**
 * @class OutputSLHAfromFH
 * @ingroup SUSY
 * @brief A class for outputting SUSY parameters to an SLHA file.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class OutputSLHAfromFH : public ThObservable {
public:

    OutputSLHAfromFH(const StandardModel& SM_i)
    : ThObservable(SM_i), output(true), mySUSY(static_cast<const SUSY*> (&SM_i))
    {
        if (mySUSY->isModelSUSY() == false)
            throw std::runtime_error("\nERROR: The SUSY mass spectrum can only be computed in a SUSY model. Please check your observables list.\n");
    };

    double computeThValue()
    {
        if (output) {
            output = false;
            (mySUSY->getMyFH()->OutputSLHA("output.slha"));
        }

        return 0.0;
    };

private:
    bool output;
    const SUSY * mySUSY;
};

#endif	/* OUTPUTSLHAFROMFH_H */

