/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef METASTABILITY_H
#define METASTABILITY_H

#include <gslpp.h>
#include <gsl/gsl_math.h>
#include "ThObservable.h"
#include "ScalarPotential.h"
#include "SUSY.h"

/*******************************************************************************
 * GSL Function Conversion BEGIN                                                  *
 * ****************************************************************************/

template<class F>
static double gslFunctionAdapter( double x, void* p)
{
    // Here I do recover the "right" pointer, safer to use static_cast
    // than reinterpret_cast.
    F* function = static_cast<F*>( p );
    return (*function)( x );
}

template<class F>
gsl_function convertToGslFunction( const F& f )
{
    gsl_function gslFunction;
    
    const void* p = &f;
    assert (p != 0);
    
    gslFunction.function = &gslFunctionAdapter<F>;
    // Just to eliminate the const.
    gslFunction.params = const_cast<void*>( p );
    
    return gslFunction;
}

/*******************************************************************************
 * GSL Function conversion END                                                     *
 * ****************************************************************************/

/**
 * @class FindAction
 * @ingroup SUSY
 * @brief FindAction.
 */
class FindAction: public ThObservable {
public:

    /**
     * @brief FindAction constructor.
     */
    FindAction(const StandardModel& SM_i);

    /**
     * @brief FindAction destructor.
     */
    ~FindAction();

    /**
     * @return @f$FindAction@f$
     */
    double computeThValue();

    double rpar, phi0par, dVpar, d2Vpar, delta_phi_cutoffpar, x1par, x2par, x3par;
    gslpp::vector<double> potentialcoefficientspar;

    gslpp::vector<double> rkqs(double y01, double y02, double dydr0, double r0, double dr, double epsfrac[2], double epsabs[2]);

protected:
    SUSYScalarPotential * mySUSYScalarPotential;

private:
    const SUSY& mySUSY;

    double func(double rpar)
    {
        double exsol=ExactSolution(rpar, phi0par, dVpar, d2Vpar)(0);
        return fabs(exsol)-fabs(delta_phi_cutoffpar);
    };

    double invertedpotential(double x)
    {
        return -mySUSYScalarPotential->potential(potentialcoefficientspar, x*x1par, x*x2par, x*x3par);
    }

    double fPS(double x);
    double Simpsonintegrand(double r, double phi, double dphi, double VphiMin_i);
    double deformedV(double phi);
    gslpp::vector<double> InitialConditions(double delta_phi0, double rmin, double delta_phi_cutoff, double distance, double dV_at_delta_phi0, double d2V_at_phi0);
//    gslpp::vector<double> splinepath(double x1, double x2, double x3, double Vmin0, gslpp::vector<double> dV0, double Vmin, gslpp::vector<double> dV);
    gslpp::vector<double> ExactSolution(double r, double phi0, double dV, double d2V);
    gslpp::vector<double> integrateProfile(double r0, double y01, double y02, double dr0, double epsfrac[2], double epsabs[2], double drmin, double rmax, double distance);
    gslpp::vector<double> dY(double y1, double y2, double r);
    int dYfunc(double r, const double y[], double ODE[], void *flags);
    int dYJac(double r, const double y[], double *dfdy, double dfdt[], void *order);

};

#endif /* METASTABILITY_H */

