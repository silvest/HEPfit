/* 
 * File:   Vud.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:45 PM
 */

#ifndef VUD_H
#define	VUD_H

#include <StandardModel.h>
#include <ThObservable.h>

class Vud : public ThObservable {
public:
    Vud(StandardModel *);
    Vud(const Vud& orig);
    virtual ~Vud();
    double getThValue();
private:
    StandardModel * myModel;
};

#endif	/* VUD_H */

