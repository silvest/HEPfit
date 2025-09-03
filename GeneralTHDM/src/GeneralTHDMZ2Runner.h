/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMZ2RUNNER_H
#define GENERALTHDMZ2RUNNER_H

#include "GeneralTHDMZ2.h"

/**
 * @class GeneralTHDMZ2Runner
 * @ingroup GeneralTHDMZ2
 * @brief An RGE running algorithm for the GeneralTHDMZ2 parameters.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Renormalization group evolution of the relevant SM and GeneralTHDMZ2 parameters
 */
class GeneralTHDMZ2Runner {
public:

    /**
     * @brief GeneralTHDMZ2Runner constructor.
     */
    GeneralTHDMZ2Runner(const StandardModel& SM_i);

    /**
     * @brief Runner destructor.
     */
    virtual ~GeneralTHDMZ2Runner();

    virtual double RGEGeneralTHDMZ2Runner(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order);

    /**
     *
     * @brief The beta function of @f$\lambda_1@f$ appearing in unitarity conditions
     * @return @f$\beta_{\lambda_1}@f$ as defined in Cacchio:2016qyh
     */
    double betalambda1_Z2(double la1, double la3, double la4, double la5, double Yb1, double Ytau1){
        return ((12.*la1*la1 + 4.*la3*la3 + 4.*la3*la4 + 2.*la4*la4 + 2.*la5*la5 +
                 12.*la1*Yb1*Yb1 + 4.*la1*Ytau1*Ytau1 - 12.*Yb1*Yb1*Yb1*Yb1 -
                 4.*Ytau1*Ytau1*Ytau1*Ytau1)/16./M_PI/M_PI);
    }

    /**
     *
     * @brief The beta function of @f$\lambda_2@f$ appearing in unitarity conditions
     * @return @f$\beta_{\lambda_2}@f$ as defined in Cacchio:2016qyh
     */
    double betalambda2_Z2(double la2, double la3, double la4, double la5, double Yt, double Yb2, double Ytau2){
        return ((12.*la2*la2 + 4.*la3*la3 + 4.*la3*la4 + 2.*la4*la4 + 2.*la5*la5 +
                 12.*la2*Yb2*Yb2 + 4.*la2*Ytau2*Ytau2 + 12.*la2*Yt*Yt - 12.*Yb2*Yb2*Yb2*Yb2 -
                 4.*Ytau2*Ytau2*Ytau2*Ytau2 - 12.*Yt*Yt*Yt*Yt)/16./M_PI/M_PI);
    }

    /**
     *
     * @brief The beta function of @f$\lambda_3@f$ appearing in unitarity conditions
     * @return @f$\beta_{\lambda_3}@f$ as defined in Cacchio:2016qyh
     */
    double betalambda3_Z2(double la1, double la2, double la3, double la4, double la5, double Yt, double Yb1, double Yb2, double Ytau1, double Ytau2){
        return ((4.*la3*la3 + 2.*la4*la4 + (la1 + la2)*(6.*la3 + 2.*la4) + 2.*la5*la5 +
                 6.*la3*Yb1*Yb1 + 6.*la3*Yb2*Yb2 + 2.*la3*Ytau1*Ytau1 + 2.*la3*Ytau2*Ytau2 +
                 6.*la3*Yt*Yt - 12.*Yb1*Yb1*(Yb2*Yb2 + Yt*Yt) - 4.*Ytau1*Ytau1*Ytau2*Ytau2)/16./M_PI/M_PI);
    }

    /**
     *
     * @brief The beta function of @f$\lambda_4@f$ appearing in unitarity conditions
     * @return @f$\beta_{\lambda_4}@f$ as defined in Cacchio:2016qyh
     */
    double betalambda4_Z2(double la1, double la2, double la3, double la4, double la5, double Yt, double Yb1, double Yb2, double Ytau1, double Ytau2){
        return (((2.*la1 + 2.*la2 + 8.*la3)*la4 + 4.*la4*la4 + 8.*la5*la5 +
                  6.*la4*Yb1*Yb1 + 6.*la4*Yb2*Yb2 + 2.*la4*Ytau1*Ytau1 + 2.*la4*Ytau2*Ytau2 +
                  6.*la4*Yt*Yt - 4.*Ytau1*Ytau1*Ytau2*Ytau2 - 12.*Yb1*Yb1*(Yb2*Yb2 - Yt*Yt))/16./M_PI/M_PI);
    }

    /**
     *
     * @brief The beta function of @f$\lambda_5@f$ appearing in unitarity conditions
     * @return @f$\beta_{\lambda_5}@f$ as defined in Cacchio:2016qyh
     */
    double betalambda5_Z2(double la1, double la2, double la3, double la4, double la5, double Yt, double Yb1, double Yb2, double Ytau1, double Ytau2){
        return (((2.*la1 + 2.*la2 + 8.*la3 + 12.*la4)*la5 +
                 6.*la5*Yb1*Yb1 + 2.*la5*Ytau1*Ytau1 + 2.*la5*(3.*Yb2*Yb2 + Ytau2*Ytau2 + 3.*Yt*Yt) -
                 12.*Yb1*Yb1*Yb2*Yb2 - 4.*Ytau1*Ytau1*Ytau2*Ytau2)/16./M_PI/M_PI);
    }

    /**
     *
     * @brief The public function which contains all relevant GTHDMZ2 parameter after running
     * @return @f$\big((\lambda_1,\lambda_2,\lambda_3,\lambda_4,\lambda_5),(Y_t,Y_{b,1},Y_{b,2},Y_{\tau,1},Y_{\tau,2}),(\text{WFR}_1,\text{WFR}_2,\text{WFR}_3,\text{WFR}_4,0)\big)@f$
     */
    gslpp::matrix<double> getGTHDMZ2_at_Q();


    const GeneralTHDMZ2 * myGTHDMZ2;

private:

    void runGeneralTHDMZ2parameters();

    void computeWFR_Z2();

    gslpp::matrix<double> myZ2_at_Q;

    double vev, cW2, Ale, Als, MZ, MZ2;
    double m12_2, mHl2, mHh2, mA2, mHp2;
    double beta, tanb, sinb, cosb, bma, alpha;
    double Q_cutoff, Q_GTHDM, Rpeps, NLOuniscale;
    double g1_at_Q, g2_at_Q, g3_at_Q;
    double Ytop_at_Q, Ybottom1_at_Q, Ybottom2_at_Q, Ytau1_at_Q, Ytau2_at_Q;
    double m11_2_at_Q, m22_2_at_Q, m12_2_at_Q;
    double lambda1_at_Q, lambda2_at_Q, lambda3_at_Q, lambda4_at_Q, lambda5_at_Q;
    double WFRcomb1, WFRcomb2, WFRcomb3, WFRcomb4;
};

#endif /* GENERALTHDMZ2RUNNER_H */

