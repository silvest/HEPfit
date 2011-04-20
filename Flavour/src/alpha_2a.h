/* 
 * File:   alpha_2a.h
 * Author: silvest
 *
 * Created on April 4, 2011, 10:58 AM
 */

#ifndef ALPHA_2A_H
#define	ALPHA_2A_H

#include <StandardModel.h>
#include <ThObservable.h>

class alpha_2a : public ThObservable {
public:
    alpha_2a(StandardModel *);
    alpha_2a(const alpha_2a& orig);
    virtual ~alpha_2a();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* ALPHA_2A_H */

