/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GG4L_H
#define	GG4L_H

#include <SusyFit.h>
#include "myModel.h"


/**
 * @class gg4l
 * @ingroup gg4l
 * @brief A class for the gg -> 4l.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class gg4l : public ThObservable {
public:
    gg4l(const StandardModel& SM_i);
    virtual ~gg4l();
    void updateParameters();
    
protected:
    
    /* Define the couplings here. */
    double ct;
    double cg;
    double cV;
    double cA;
    double sw2;
    double fact;
    double kfact;

private:
    const myModel * my_model;

};

class yield : public gg4l {
public:
    
    /**
     * @brief yield
     */
    yield(const StandardModel& SM_i, unsigned int bin_i);
    
    /**
     * @return return the BR
     */
    double computeThValue ();
    
private:
    unsigned int bin;

};

class C_V : public gg4l{
public:
    
    /**
     * @brief \f$ C_V \f$
     */
    C_V(const StandardModel& SM_i);
    
    /**
     * @return return the BR
     */
    double computeThValue ();
    
private:
    
};

class C_A : public gg4l{
public:
    
    /**
     * @brief \f$ C_A \f$
     */
    C_A(const StandardModel& SM_i);
    
    /**
     * @return return the BR
     */
    double computeThValue ();
    
private:
    
};

#endif	/* GG4LL_H */

