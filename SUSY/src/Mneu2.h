/* 
 * File:   Mneu2.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.33
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

