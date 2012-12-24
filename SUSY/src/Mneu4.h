/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MNEU4_H
#define	MNEU4_H

class Mneu4 : public ThObservable {
public:
    Mneu4(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return ((static_cast<const SUSY*> (&SM))->getMneu()(3));
    };
private:

};

#endif	/* MNEU4_H */

