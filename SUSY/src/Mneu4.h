/* 
 * File:   Mneu4.h
 * Author: girardi_mac
 *
 * Created on 28 settembre 2012, 16.33
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

