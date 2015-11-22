/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_TEST_H
#define	EW_TEST_H

#include "ThObservable.h"
#include "StandardModel.h"
#include "EW_ABC.h"
#include "EW_BURGESS.h"
#include "EW_CHMN.h"
#include <string>

/**
 * @class EW_TEST
 * @ingroup NewPhysics
 * @brief A test class for the EWPO.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EW_TEST : public ThObservable {
public:

    EW_TEST(const std::string mode_i, const std::string obsname_i, const StandardModel& SM_i);

    double computeThValue();

private:
    const std::string mode;
    const std::string obsname;
    const StandardModel& SM;
    const EW_ABC *myEW_ABC;
    const EW_BURGESS *myEW_BURGESS;
    const EW_CHMN *myEW_CHMN;

};

#endif	/* EW_TEST_H */

