/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
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
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class OutputSLHAfromFH : public ThObservable {
public:

    OutputSLHAfromFH(const StandardModel& SM_i)
    : ThObservable(SM_i), output(true)
    {
    };

    double computeThValue()
    {
        if (output) {
            output = false;
            ((static_cast<const SUSY*> (&SM))->getMyFH()->OutputSLHA("output.slha"));
        }

        return 0.0;
    };

private:
    bool output;

};

#endif	/* OUTPUTSLHAFROMFH_H */

