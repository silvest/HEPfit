/* 
 * File:   Flavour.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:49 PM
 */

#ifndef FLAVOUR_H
#define	FLAVOUR_H

#include <ThObsType.h>
#include <StandardModel.h>
#include <DmBd.h>
#include <DmBd0.h>
#include <Vud.h>
#include <Vus.h>
#include <Vcb.h>
#include <Vub.h>
#include <alpha.h>
#include <alpha_2a.h>
#include <gamma.h>

using namespace gslpp;

class Flavour : public ThObsType {
public:

    Flavour(const StandardModel& SM_i) : ThObsType(SM_i), DF2(SM_i) {
    };

    Flavour(const Flavour& orig) : ThObsType(orig.SM), DF2(orig.SM) {
     };

    virtual ~Flavour() {
    };

private:
    HeffDF2 DF2;
};

#endif	/* FLAVOUR_H */
