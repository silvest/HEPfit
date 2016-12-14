/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMUNITARITY_H
#define GENERALTHDMUNITARITY_H

#include "ThObservable.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"
#include <gslpp.h>

/**
 * 
 */
class unitarity_GTHDM {
public:
    
    /**
     * @brief Constructor.
     */
    unitarity_GTHDM(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~unitarity_GTHDM();
    
    /**
     * @brief Computes the eigenvalues of the S matrix with @f$Y=2,\sigma=1@f$
     */
    bool CalcSeigen21(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i);
    
    /**
     * @brief Computes the eigenvalues of the S matrix with @f$Y=0,\sigma=1@f$
     */
    bool CalcSeigen01(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i);
    
    /**
     * @brief Computes the eigenvalues of the S matrix with @f$Y=0,\sigma=0@f$
     */
    bool CalcSeigen00(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i);
    
    /**
     * @brief Assigns to a vector the eigenvalues of the S matrix with @f$Y=2,\sigma=1@f$
     */
    gslpp::vector<double> getSeigen21();
    
    /**
     * @brief Assigns to a vector the eigenvalues of the S matrix with @f$Y=0,\sigma=1@f$
     */
    gslpp::vector<double> getSeigen01();
    
    /**
     * @brief Assigns to a vector the eigenvalues of the S matrix with @f$Y=0,\sigma=0@f$
     */
    gslpp::vector<double> getSeigen00();
    
    
private:
    const GeneralTHDM& myGTHDM;
    
    gslpp::matrix<gslpp::complex> Smat21, Smat01, Smat00; 
    gslpp::matrix<gslpp::complex> Seigvec21, Seigvec01, Seigvec00;
    gslpp::vector<double> Seigval21, Seigval01, Seigval00;
};


/**
 * 
 */
class unitarity1_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity1_GTHDM constructor.
     */
    unitarity1_GTHDM(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};


/**
 * 
 */
class unitarity2_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity2_GTHDM constructor.
     */
    unitarity2_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity3_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity3_GTHDM constructor.
     */
    unitarity3_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity5_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity5_GTHDM constructor.
     */
    unitarity5_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity6_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity6_GTHDM constructor.
     */
    unitarity6_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity7_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity7_GTHDM constructor.
     */
    unitarity7_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity8_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity8_GTHDM constructor.
     */
    unitarity8_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity9_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity9_GTHDM constructor.
     */
    unitarity9_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity10_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity10_GTHDM constructor.
     */
    unitarity10_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity11_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity11_GTHDM constructor.
     */
    unitarity11_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity12_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity12_GTHDM constructor.
     */
    unitarity12_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    unitarity_GTHDM myunitarity_GTHDM;
};

/**
 * 
 */
class unitarity4_GTHDM: public ThObservable {
public:

    /**
     * @brief unitarity4_GTHDM constructor.
     */
    unitarity4_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

#endif /* GENERALTHDMUNITARITY_H */

