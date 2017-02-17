/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LRSMQUANTITIES_H
#define	LRSMQUANTITIES_H

#include <stdexcept>
#include "ThObservable.h"
#include "LeftRightSymmetricModel.h"

//class LeftRightSymmetricModel;

/**
 * @class LRSMquantities
 * @ingroup LeftRightSymmetricModel
 * @brief A class for calculating the Higgs spectrum and other potential parameters at tree level.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class LRSMquantities {

public:

    LRSMquantities(const LeftRightSymmetricModel & LRSM_in);

    /**
     * @brief LRSMquantities constructor.
     * @details 
     */
//    LRSMquantities(const StandardModel& SM_i);

    /**
     * @brief LRSMquantities destructor.
     * @details 
     */
    ~LRSMquantities();
    

    /**
     * @brief Computes the exact neutral spectrum at tree level.
     */
    bool CalcNeutralMasses(gslpp::matrix<gslpp::complex>& U_i, double mH0sq[5]);

    /**
     * @brief Computes the approximate neutral spectrum at tree level.
     */
    bool CalcNeutralMasses_app(double mH0sq_app[4]);

private:
    const LeftRightSymmetricModel& myLRSM;
//    const LeftRightSymmetricModel * myLRSM;

    /**
     * @brief Stores the tree-level neutral mass square matrix.
     */
    gslpp::matrix<double> Msqneutral;
    
    /**
     * @brief Stores the tree-level neutral mixing matrix.
     */
//    gslpp::matrix<gslpp::complex> U;

    /**
     * @brief Stores the tree-level mass square eigenvalues.
     */
//    gslpp::matrix<gslpp::complex> mH0sq;
};

/**
 * @class mu1_2_LRSM
 * @ingroup LRSM
 * @brief 
 */
class mu1_2_LRSM: public ThObservable {
public:

    /**
     * @brief mu1_2_LRSM constructor.
     */
    mu1_2_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\mu_1^2@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class mu2_2_LRSM
 * @ingroup LRSM
 * @brief 
 */
class mu2_2_LRSM: public ThObservable {
public:

    /**
     * @brief mu2_2_LRSM constructor.
     */
    mu2_2_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\mu_2^2@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class mu3_2_LRSM
 * @ingroup LRSM
 * @brief 
 */
class mu3_2_LRSM: public ThObservable {
public:

    /**
     * @brief mu3_2_LRSM constructor.
     */
    mu3_2_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\mu_3^2@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class rho2_LRSM
 * @ingroup LRSM
 * @brief 
 */
class rho2_LRSM: public ThObservable {
public:

    /**
     * @brief rho2_LRSM constructor.
     */
    rho2_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\rho_2@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class rho3_LRSM
 * @ingroup LRSM
 * @brief 
 */
class rho3_LRSM: public ThObservable {
public:

    /**
     * @brief rho3_LRSM constructor.
     */
    rho3_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\rho_3@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class alpha3_LRSM
 * @ingroup LRSM
 * @brief 
 */
class alpha3_LRSM: public ThObservable {
public:

    /**
     * @brief alpha3_LRSM constructor.
     */
    alpha3_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$\alpha_3@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class MH05_LRSM
 * @ingroup LRSM
 * @brief 
 */
class MH05_LRSM: public ThObservable {
public:

    /**
     * @brief MH05_LRSM constructor.
     */
    MH05_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^0_5}@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class MH06_LRSM
 * @ingroup LRSM
 * @brief 
 */
class MH06_LRSM: public ThObservable {
public:

    /**
     * @brief MH06_LRSM constructor.
     */
    MH06_LRSM(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^0_6}@f$
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class MH01_app1
 * @ingroup LRSM
 * @brief 
 */
class MH01_app1: public ThObservable {
public:

    /**
     * @brief MH01_app1 constructor.
     */
    MH01_app1(const StandardModel& SM_i);

    /**
     * @return approximate @f$m_{H^0_1}@f$ according to 1612.09146
     */
    double computeThValue();

    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class MH0_LRSM
 * @ingroup LeftRightSymmetricModel
 * @brief A class for the scalar masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class MH0_LRSM : public ThObservable {
public:

    MH0_LRSM(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
    {
    };

    double computeThValue()
    {
        switch(index) {
            case 0:
                return sqrt(myLRSM->getmH0sq1());
            case 1:
                return sqrt(myLRSM->getmH0sq2());
            case 2:
                return sqrt(myLRSM->getmH0sq3());
            case 3:
                return sqrt(myLRSM->getmH0sq4());
            case 4:
                return sqrt(myLRSM->getmH0sq5());
            default:
                throw std::runtime_error("MH0_LRSM::computeThValue(): undefined index");
        }
    };

private:
    const int index;
    const LeftRightSymmetricModel * myLRSM;
};

/**
 * @class MH0_app
 * @ingroup LeftRightSymmetricModel
 * @brief A class for the approximate scalar masses.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class MH0_app : public ThObservable {
public:

    MH0_app(const StandardModel& SM_i, const int ind)
    : ThObservable(SM_i), index(ind), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
    {
    };

    double computeThValue()
    {
        switch(index) {
            case 0:
                return sqrt(myLRSM->getmH0sq1_app());
            case 1:
                return sqrt(myLRSM->getmH0sq2_app());
            case 2:
                return sqrt(myLRSM->getmH0sq3_app());
            case 3:
                return sqrt(myLRSM->getmH0sq4_app());
            default:
                throw std::runtime_error("MH0_app::computeThValue(): undefined index");
        }
    };

private:
    const int index;
    const LeftRightSymmetricModel * myLRSM;
};

#endif	/* LRSMQUANTITIES_H */
