/* 
 * File:   Dmb.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:51 PM
 */

#ifndef DMB_H
#define	DMB_H

#include <StandardModel.h>
#include <ThObservable.h>
#include <gslpp_complex.h>

class Dmb : public ThObservable {
public:
    Dmb(StandardModel *, const int);
    Dmb(const Dmb& orig);
    virtual ~Dmb();
    double getThValue();
private:
    StandardModel * myModel;
    int LE;
};

#endif	/* DMB_H */

