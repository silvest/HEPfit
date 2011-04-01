/* 
 * File:   Flavour.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:49 PM
 */

#ifndef FLAVOUR_H
#define	FLAVOUR_H

#include "Dmb.h"
#include "Vud.h"
#include "Vus.h"
#include "Vcb.h"
#include "Vub.h"
#include "alpha.h"
#include "gamma.h"

class Flavour {
public:
    Flavour(StandardModel&);
    Flavour(const Flavour& orig);
    virtual ~Flavour();
private:
    StandardModel * myModel;
};

#endif	/* FLAVOUR_H */

