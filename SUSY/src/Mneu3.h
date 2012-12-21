/* 
 * File:   Mneu3.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.33
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

