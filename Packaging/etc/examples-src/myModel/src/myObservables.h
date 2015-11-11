/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYOBSERVABLES_H
#define	MYOBSERVABLES_H

#include <HEPfit.h>
#include "myModel.h"


/**
 * @class myObservables
 * @brief A class for the gg -> 4l.
 */
class myObservables : public ThObservable {
public:
    myObservables(const StandardModel& SM_i);
    virtual ~myObservables();
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

class yield : public myObservables {
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

class C_3 : public myObservables {
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

class C_4 : public myObservables {
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

#endif	/* MYOBSERVABLES_H */

