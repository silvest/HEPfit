/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THOBSTYPE_H
#define	THOBSTYPE_H

#include <StandardModel.h>

class ThObsType {
public:
    ThObsType(const StandardModel& SM_i) : SM(SM_i) {};
    ThObsType(const ThObsType& orig) : SM(orig.SM) {};
    virtual ~ThObsType() {};
    const StandardModel& getModel() const {
        return SM;
    };

protected:
    const StandardModel& SM;
};

#endif	/* THOBSTYPE_H */
