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
 * @class Rlow_BXqll
 * @ingroup Flavour
 * @brief A class for @f$R_{low}(\hat s)@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute @f$Rlow_{quark}(\hat s)@f$ as a function
 * of @f$\hat s = \frac{q^2}{m_b^2}@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * Rlow_{quark}(\hat s)= \int{0.05}{0.25} \frac{1}{\Gamma(B\to X_c e\nu}\frac {d\Gamma}{d\hat s} \,.
 * @f]
 */
class Rlow_BXqll : public ThObservable{
public:
    
     /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Rlow_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
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
 * @class Rhigh_BXqll
 * @ingroup Flavour
 * @brief A class for @f$R_{high}(\hat s)@f$ in @f$B \to X_q l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute @f$Rhigh_{quark}(\hat s)@f$ as a function
 * of @f$\hat s = \frac{q^2}{m_b^2}@f$ for @f$B \to X_q l^+l^-@f$ :
 * @f[
 * Rhigh_{quark}(\hat s)= \int{0.6}{1.0} \frac{1}{\Gamma(B\to X_c e\nu}\frac {d\Gamma}{d\hat s} \,.
 * @f]
 */
class Rhigh_BXqll : public ThObservable{
public:
    
     /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    Rhigh_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<BR>@f$ in @f$M \to P l^+l^-@f$.
    * @return @f$<BR>@f$
    */
    double computeThValue ();
    
private:
    QCD::lepton lep; /**< Final leptons type. */
    QCD::quark quark; /**< Final quark type. */
//    BXqll myBXqll;
};



#endif	/* BXQLLOBSERVABLES_H */
