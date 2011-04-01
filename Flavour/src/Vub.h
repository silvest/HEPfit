/* 
 * File:   Vub.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:41 PM
 */

#ifndef VUB_H
#define	VUB_H

#include <StandardModel.h>
#include <ThObservable.h>

class Vub : public ThObservable {
public:
    Vub(StandardModel *);
    Vub(const Vub& orig);
    virtual ~Vub();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* VUB_H */

