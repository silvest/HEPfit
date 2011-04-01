/* 
 * File:   alpha.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef ALPHA_H
#define	ALPHA_H

#include <StandardModel.h>
#include <ThObservable.h>

class alpha : public ThObservable {
public:
    alpha(StandardModel *);
    alpha(const alpha& orig);
    virtual ~alpha();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* ALPHA_H */

