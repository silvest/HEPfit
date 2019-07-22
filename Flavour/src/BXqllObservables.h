/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BXQLLOBSERVABLES_H
#define	BXQLLOBSERVABLES_H

#include "BXqll.h"
#include "ThObservable.h"

/**
 * @class R_BXqll
 * @ingroup Flavour
 * @brief A class for @f$R(\hat s)@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute @f$R_{quark}(\hat s)@f$ as a function
 * of @f$\hat s = \frac{q^2}{m_b^2}@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * R_{quark}(\hat s)= \frac{1}{\Gamma(B\to X_c e\nu}\frac {d\Gamma}{d\hat s} \,.
 * @f]
 */
class R_BXqll : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    R_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<BR>@f$ in @f$M \to P l^+l^-@f$.
    * @return @f$<BR>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::quark quark; /**< Final quark type. */
    BXqll myBXqll;
};


/**
 * @class HT_BXqll
 * @ingroup Flavour
 * @brief A class for @f$H_T@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integral of @f$H_T(q^2)@f$
 * in bins of @f$q^2@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * H_T[q_m^2,q_M^2] = \int{q_m^2}{q_M^2} dq^2 \int{-1}{+1} dz \frac{d^2\Gamma}{dq^2 dz^2}
 * \frac{2}{3} P_0(z) + \frac{10}{3} P_2(z) \,.
 * @f]
 */
class HT_BXqll : public ThObservable{
public:
    
     /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    HT_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$H_T@f$ in @f$B \to X_q l^+l^-@f$.
    * @return @f$H_T@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::quark quark; /**< Final quark type. */
    BXqll myBXqll;
};


/**
 * @class HL_BXqll
 * @ingroup Flavour
 * @brief A class for @f$H_L@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integral of @f$H_L(q^2)@f$
 * in bins of @f$q^2@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * H_L[q_m^2,q_M^2] = \int{q_m^2}{q_M^2} dq^2 \int{-1}{+1} dz \frac{d^2\Gamma}{dq^2 dz^2}
 * \frac{1}{3} P_0(z) - \frac{10}{3} P_2(z) \,.
 * @f]
 */
class HL_BXqll : public ThObservable{
public:
    
     /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    HL_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$H_L@f$ in @f$B \to X_q l^+l^-@f$.
    * @return @f$H_L@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::quark quark; /**< Final quark type. */
    BXqll myBXqll;
};


/**
 * @class HA_BXqll
 * @ingroup Flavour
 * @brief A class for @f$H_A@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the integral of @f$H_A(q^2)@f$
 * in bins of @f$q^2@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * H_A[q_m^2,q_M^2] = \int{q_m^2}{q_M^2} dq^2 \int{-1}{+1} dz \frac{d^2\Gamma}{dq^2 dz^2}
 * \frac{4}{3} {\rm sign}(z) \,.
 * @f]
 */
class HA_BXqll : public ThObservable{
public:
    
     /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    HA_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$H_A@f$ in @f$B \to X_q l^+l^-@f$.
    * @return @f$H_A@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::quark quark; /**< Final quark type. */
    BXqll myBXqll;
};


#endif	/* BXQLLOBSERVABLES_H */
