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
    double c1;
    double c2;
    double c3;
    double c4;
    double sw2;
    double fact;
    double kfact;

private:
    const myModel * my_model;

};

class yield : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    yield(const StandardModel& SM_i, unsigned int bin_i);
    
    /**
     * @return yield
     */
    double computeThValue ();
    
private:
    unsigned int bin;

};

class C_3 : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    C_3(const StandardModel& SM_i);
    
    /**
     * @return C_V
     */
    double computeThValue ();
    
private:
    
};

class C_4 : public gg4l {
public:
    
    /**
     * @brief Constructor.
     */
    C_4(const StandardModel& SM_i);
    
    /**
     * @return C_A
     */
    double computeThValue ();
    
private:
    
};

#endif	/* GG4LL_H */

