/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELOBSERVABLE_H
#define	MODELOBSERVABLE_H

#include "ThObsType.h"
#include "StandardModel.h"
#include <stdexcept>

/**
 * @class ModelObservable
 * @ingroup Observable
 * @brief A class for model parameters. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class ModelObservable : public ThObsType {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel()
     */
    ModelObservable(const StandardModel & SM_i)
    : ThObsType(SM_i)
    {};
    
    /**
     * @brief The default destructor.
     */
    virtual ~ModelObservable()
    {};
};

#endif	/* MODELOBSERVABLE_H */

