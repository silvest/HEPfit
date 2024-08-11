/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEPTON_DEC_H
#define	LEPTON_DEC_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */

/**
 * @class GammaMuon 
 * @ingroup EW
 * @brief An observable class for the muon decay width
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the muon decay width 
 * @f$\Gamma(\mu \to e \nu \bar{\nu}@f$
 *
 */
class GammaMuon : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaMuon(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The muon decay width @f$\Gamma(\mu \to e \nu \bar{\nu}@f$
     * @return @f$\Gamma(\mu \to e \nu \bar{\nu}@f$
     */
    double computeThValue();

private:

};


/**
 * @class GammaTautoMuon 
 * @ingroup EW
 * @brief An observable class for the leptonic tau decay width into muons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau decay width 
 * @f$\Gamma(\tau \to \mu \nu \bar{\nu})@f$
 *
 */
class GammaTautoMuon : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaTautoMuon(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The leptonic tau decay width @f$\Gamma(\tau \to \mu \nu \bar{\nu})@f$
     * @return @f$\Gamma(\tau \to \mu \nu \bar{\nu})@f$
     */
    double computeThValue();

private:

};


/**
 * @class GammaTautoElectron
 * @ingroup EW
 * @brief An observable class for the leptonic tau decay width into electrons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau decay width 
 * @f$\Gamma(\tau \to e \nu \bar{\nu})@f$
 *
 */
class GammaTautoElectron : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaTautoElectron(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The leptonic tau decay width @f$\Gamma(\tau \to e \nu \bar{\nu})@f$
     * @return @f$\Gamma(\tau \to e \nu \bar{\nu})@f$
     */
    double computeThValue();

private:

};


// LFU Tests in Tau decays


/**
 * @class gmugeTauLFU
 * @ingroup EW
 * @brief An observable class for the tau LFU ratio @f$g_\mu/g_e@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau LFU ratio obtained from
 * @f$\frac{\tau \to \mu \nu \nu }{Gamma(\tau \to e \nu \nu } @f$
 *
 */
class gmugeTauLFU : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gmugeTauLFU(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The tau LFU ratio @f$g_\mu/g_e@f$
     * @return @f$g_\mu/g_e@f$
     */
    double computeThValue();

private:

};


/**
 * @class gtaugmuTauLFU
 * @ingroup EW
 * @brief An observable class for the tau LFU ratio @f$g_\tau/g_\mu@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau LFU ratio obtained from
 * @f$\frac{\tau \to e \nu \nu}{\mu \to e \nu \nu } @f$
 *
 */
class gtaugmuTauLFU : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gtaugmuTauLFU(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The tau LFU ratio @f$g_\tau/g_\mu@f$
     * @return @f$g_\tau/g_\mu@f$
     */
    double computeThValue();

private:

};


/**
 * @class gtaugeTauLFU
 * @ingroup EW
 * @brief An observable class for the tau LFU ratio @f$g_\tau/g_e@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau LFU ratio obtained from
 * @f$\frac{\tau \to \mu \nu \nu }{ \mu \to e \nu \nu } @f$
 *
 */
class gtaugeTauLFU : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gtaugeTauLFU(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The tau LFU ratio @f$g_\tau/g_e@f$
     * @return @f$g_\tau/g_e@f$
     */
    double computeThValue();

private:

};


/**
 * @class gtaugmuPiTauLFU
 * @ingroup EW
 * @brief An observable class for the tau LFU ratio @f$\left(g_\tau/ g_\mu\right)_\pi @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau LFU ratio obtained from
 * @f$\frac{\tau \to \pi \nu}{\pi \to \mu \nu } @f$
 *
 */
class gtaugmuPiTauLFU : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gtaugmuPiTauLFU(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The tau LFU ratio @f$\left(g_\tau/ g_\mu\right)_\pi @f$
     * @return @f$\left(g_\tau/ g_\mu\right)_\pi @f$
     */
    double computeThValue();

private:

};


/**
 * @class gtaugmuKTauLFU
 * @ingroup EW
 * @brief An observable class for the tau LFU ratio @f$\left(g_\tau/ g_\mu\right)_K @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the tau LFU ratio obtained from
 * @f$\frac{\tau \to K \nu}{K \to \mu \nu } @f$
 *
 */
class gtaugmuKTauLFU : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gtaugmuKTauLFU(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The tau LFU ratio @f$\left(g_\tau/ g_\mu\right)_K @f$
     * @return @f$\left(g_\tau/ g_\mu\right)_K @f$
     */
    double computeThValue();

private:

};




/** 
 * @}
 */

#endif	/* LEPTON_DEC_H */
