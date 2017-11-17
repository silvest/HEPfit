/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PMNS_H
#define	PMNS_H

#include <math.h>
#include <gslpp.h>

/**
 * @class PMNS
 * @ingroup StandardModel
 * @brief A class for the %PMNS matrix elements.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining the %PMNS matrix and its elements with
 * the Wolfenstein and Gilman parameterizations. 
 */
class PMNS {
public:
    PMNS();
    ~PMNS();

    void setPMNS(double s12_v, double s13_v, double s23_v, double delta_v, double alpha21_v, double alpha31_v);

    gslpp::matrix<gslpp::complex> getPMNS() const {
        return U;
    }

    void setPMNS(gslpp::matrix<gslpp::complex> U) {
        this->U = U;
    }

    // Gilman parameterization
    double gets12();
    double gets13();
    double gets23();
    double getc12();
    double getc13();
    double getc23();
    double getdelta();
    double getalpha21();
    double getalpha31();

    //Absolute values of PMNS elements
    double getUe1();
    double getUe2();
    double getUe3();
    double getUmu1();
    double getUmu2();
    double getUmu3();
    double getUtau1();
    double getUtau2();
    double getUtau3();

    //Phases of PMNS elements
    double getArgUe1();
    double getArgUe2();
    double getArgUe3();
    double getArgUmu1();
    double getArgUmu2();
    double getArgUmu3();
    double getArgUtau1();
    double getArgUtau2();
    double getArgUtau3();

    //Complex values of PMNS elements
    gslpp::complex U_e1();
    gslpp::complex U_e2();
    gslpp::complex U_e3();
    gslpp::complex U_mu1();
    gslpp::complex U_mu2();
    gslpp::complex U_mu3();
    gslpp::complex U_tau1();
    gslpp::complex U_tau2();
    gslpp::complex U_tau3();

    
private:
    double s12, s13, s23, delta, alpha21, alpha31;
    double c12, c23, c13;

    gslpp::matrix<gslpp::complex> U;

};

#endif	/* PMNS_H */

