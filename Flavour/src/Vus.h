/* 
 * File:   Vus.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:17 PM
 */

#ifndef VUS_H
#define	VUS_H

#include <StandardModel.h>
#include <ThObservable.h>

class Vus : public ThObservable {
public:
    Vus(StandardModel *);
    Vus(const Vus& orig);
    virtual ~Vus();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* VUS_H */

