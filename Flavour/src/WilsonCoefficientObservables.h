/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONCOEFFICIENTOBSERVABLES_H
#define WILSONCOEFFICIENTOBSERVABLES_H

class StandardModel;
#include "ThObservable.h"
#include "QCD.h"

/**
 * @class WC_C7g
 * @ingroup Flavour
 * @brief A class for the Wilson coefficient @f$C_7@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient @f$C_7(\mu)@f$ 
 * at a specific scale @f$\mu@f$.
 */
class WC_C7g : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] part_i an unsigned integer used as toggle for real, imaginary, absolute or argument value.
     */
    WC_C7g(const StandardModel& SM_i, unsigned int part_i);
    
    /**
    * @brief The Wilson coefficient @f$C_7(\mu)@f$.
    * @return @f$C_7(\mu) @f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    double mu; /**< The scale for the Wilson coefficient. */
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    unsigned int part;/**< toggle for real, imaginary, absolute or argument value */
};

/**
 * @class WC_C9
 * @ingroup Flavour
 * @brief A class for the Wilson coefficient @f$C_9^{lep}(\mu)@f$ for a specific leptonic final state. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient @f$C_9^{lep}(\mu)@f$ 
 * at a specific scale @f$\mu@f$ for a specific leptonic final state.
 */
class WC_C9 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] part_i an unsigned integer used as toggle for real, imaginary, absolute or argument value.
     * @param[in] lep_i final leptons of the decay
     */
    WC_C9(const StandardModel& SM_i, unsigned int part_i, QCD::lepton lep_i);
    
    /**
    * @brief The Wilson coefficient @f$C_9^{lep}(\mu)@f$.
    * @return @f$C_9^{lep}(\mu)@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lepton; /**< Final leptons type. */
    const StandardModel& mySM;
    double mu; /**< The scale for the Wilson coefficient. */
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    unsigned int part;/**< toggle for real, imaginary, absolute or argument value */
};

/**
 * @class WC_C10
 * @ingroup Flavour
 * @brief A class for the Wilson coefficient @f$C_10^{lep}(\mu)@f$ for a specific leptonic final state. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient @f$C_10^{lep}(\mu)@f$ 
 * at a specific scale @f$\mu@f$ for a specific leptonic final state.
 */
class WC_C10 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] part_i an unsigned integer used as toggle for real, imaginary, absolute or argument value.
     * @param[in] lep_i final leptons of the decay
     */
    WC_C10(const StandardModel& SM_i, unsigned int part_i, QCD::lepton lep_i);
    
    /**
    * @brief The Wilson coefficient @f$C_10^{lep}(\mu)@f$.
    * @return @f$C_10^{lep}(\mu)@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lepton; /**< Final leptons type. */
    const StandardModel& mySM;
    double mu; /**< The scale for the Wilson coefficient. */
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    unsigned int part;/**< toggle for real, imaginary, absolute or argument value */
};

class WC_epspOeps : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] part_i an unsigned integer used as toggle for real, imaginary, absolute or argument value.
     * @param[in] mu_z scale for the Wilson coeff
     */
    WC_epspOeps(const StandardModel& SM_i, unsigned int part_i, double mu_z);
    
    /**
    * @brief The Wilson coefficients of eps' / eps.
    * @return @f$WC_epspOeps@f$
    */
    double computeThValue ();
    
private:
    const StandardModel& mySM;
    double mu; /**< The scale for the Wilson coefficient. */
    gslpp::vector<gslpp::complex> ** allcoeffv;/**<vector that contains the Wilson coeffients v */
    gslpp::vector<gslpp::complex> ** allcoeffz;/**<vector that contains the Wilson coeffients z */
    unsigned int part;/**< toggle for real, imaginary, absolute or argument value */
};

#endif	/* WILSONCOEFFICIENTOBSERVABLES_H */

    

