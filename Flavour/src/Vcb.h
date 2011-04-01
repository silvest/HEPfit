/* 
 * File:   Vcb.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef VCB_H
#define	VCB_H

#include <StandardModel.h>
#include <ThObservable.h>

class Vcb : public ThObservable {
public:
    Vcb(StandardModel *);
    Vcb(const Vcb& orig);
    virtual ~Vcb();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* VCB_H */

