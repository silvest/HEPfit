/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MNEU3_H
#define	MNEU3_H

class Mneu3 : public ThObservable {
public:
    Mneu3(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return ((static_cast<const SUSY*> (&SM))->getMneu()(2));
    };
private:

};

#endif	/* MNEU3_H */

