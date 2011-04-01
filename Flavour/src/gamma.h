/* 
 * File:   gamma.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:46 PM
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include <StandardModel.h>
#include <ThObservable.h>

class gammac : public ThObservable {
public:
    gammac(StandardModel *);
    gammac(const gammac& orig);
    virtual ~gammac();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* GAMMA_H */

