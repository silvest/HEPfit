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

class ModelObservable : public ThObsType {
public:
    ModelObservable(const StandardModel & SM_i): ThObsType(SM_i){};
};

#endif	/* MODELOBSERVABLE_H */

