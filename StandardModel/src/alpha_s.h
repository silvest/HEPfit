/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHA_S_H
#define ALPHA_S_H

#include "ThObservable.h"
#include "OrderScheme.h"
class StandardModel;

class alpha_s : public ThObservable {
public:
    /**
     * constructor
     * @param SM_i a reference to Standard Model
     * @param order the order of the computation (LO, NLO, FULLNLO, NNLO or FULLNNLO)
     */
    alpha_s(const StandardModel& SM_i, orders order);
    
    /**
     * default destructor
     */
    virtual ~alpha_s() {};
    
    /**
     * 
     * @brief A method to get the value of \f$ \alpha_{s}(\mu, order) \f$
     * @return theoretical value of \f$ \alpha_{s}(\mu, order) \f$
     */
    double computeThValue();
private:
    
    orders order;

};

#endif /* ALPHA_S_H */