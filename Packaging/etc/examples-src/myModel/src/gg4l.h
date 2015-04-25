/* 
 * Copyright (C) 2015 SusyFit Collaboration
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
 * @brief A class for the gg -> 4l.
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

class yeild : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    yeild(const StandardModel& SM_i, unsigned int bin_i);
    
    /**
     * @return yeild
     */
    double computeThValue ();
    
private:
    unsigned int bin;

};

class C_V : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    C_V(const StandardModel& SM_i);
    
    /**
     * @return C_V
     */
    double computeThValue ();
    
private:
    
};

class C_A : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    C_A(const StandardModel& SM_i);
    
    /**
     * @return C_A
     */
    double computeThValue ();
    
private:
    
};

#endif	/* GG4LL_H */

