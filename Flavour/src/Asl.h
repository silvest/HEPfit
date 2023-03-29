/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ASL_H
#define ASL_H
#include "StandardModel.h"
#include "ThObservable.h"
#include "AmpDB2.h"

#include "gslpp.h"


class Asl : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] lep_i final leptons of the decay
     */
    Asl(const StandardModel& SM_i, QCD::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl();
    
    double computeThValue ();

    private:
    QCD::lepton lep;
};
 
#endif /* ASL_H */

