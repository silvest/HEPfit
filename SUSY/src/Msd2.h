/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MD2_H
#define	MD2_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "SUSY.h"

class Msd2 : public ThObservable {
public:
    Msd2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return (sqrt((static_cast<const SUSY*> (&SM))->getMsd2()(1)));
    };
private:

};

#endif	/* MD2_H */

