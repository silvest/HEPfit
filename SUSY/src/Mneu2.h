/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MNEU2_H
#define	MNEU2_H

class Mneu2 : public ThObservable {
public:
    Mneu2(const ThObsType& ObsType) : ThObservable(ObsType) {};
    double getThValue(){
        
        return ((static_cast<const SUSY*> (&SM))->getMneu()(1));
    };
private:

};

#endif	/* MNEU2_H */

