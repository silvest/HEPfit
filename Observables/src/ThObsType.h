/* 
 * File:   ThObsType.h
 * Author: marco
 *
 * Created on April 29, 2011, 11:43 AM
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
