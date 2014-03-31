/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MUGGH_H
#define	MUGGH_H
#include <ThObservable.h>

/**
 * @class muggH
 * @ingroup HiggsCouplings
 * @brief A class for computing the ratio @f[\mu_{ggH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ggH}@f] between the gluon-gluon fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */

class muggH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muggH(const StandardModel& SM_i) : SM(SM_i) 
    {
        
    };
    muggH(const muggH& orig);
    virtual ~muggH();
    
    /**
     * method to compute the value of  @f[\mu_{ggH}@f] in the current model
     * @return 
     */
    double computeThValue();
    
private:
    StandardModel& SM;
};

#endif	/* MUGGH_H */

