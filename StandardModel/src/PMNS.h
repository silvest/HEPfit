/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PMNS_H
#define	PMNS_H

#include <math.h>
#include <gslpp.h>

using namespace gslpp;

/**
 * @class PMNS
 * @ingroup StandardModel
 * @brief A class for the PMNS matrix elements.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class PMNS {
public:
    PMNS(std::string H_i);
    PMNS(const PMNS&);
    ~PMNS();

    void setPMNS(double, double, double, double);

    void getPMNS(matrix<complex> &) const;

    // Gilman parameterization
    double gets12();
    double gets13();
    double gets23();
    double getc12();
    double getc13();
    double getc23();
    double getdelta();

    std::string getHierarchy();

    // J_CP
    double getJcp();

    // Sides
    /*double getRt(std::string H);
    double getRts(std::string H);
    double getRb(std::string H);*/

private:
    double s12;
    double s13;
    double s23;
    double delta;
    double c12;
    double c23;
    double c13;

    complex U11;
    complex U12;
    complex U13;
    complex U21;
    complex U22;
    complex U23;
    complex U31;
    complex U32;
    complex U33;

    std::string H;

};



#endif	/* PMNS_H */

