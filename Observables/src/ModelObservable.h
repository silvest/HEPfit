/* 
 * File:   ModelObservable.h
 * Author: silvest
 *
 * Created on September 28, 2012, 1:57 PM
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

