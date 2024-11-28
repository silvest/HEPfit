/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef H_DECAYS_H
#define	H_DECAYS_H

#include <stdexcept>
#include <ThObservable.h>


/**
 * @class Htobb 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to bb
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to bb 
 * @f$\Gamma(H\to bb)@f$ in GeV
 *
 */
class Htobb: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htobb(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to bb @f$\Gamma(H\to bb)@f$
     * @return @f$\Gamma(H\to bb)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htocc 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to cc
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to cc
 * @f$\Gamma(H\to cc)@f$ in GeV
 *
 */
class Htocc: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htocc(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to cc @f$\Gamma(H\to cc)@f$
     * @return @f$\Gamma(H\to cc)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htoss 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to ss
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to ss
 * @f$\Gamma(H\to ss)@f$ in GeV
 *
 */
class Htoss: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htoss(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to ss @f$\Gamma(H\to ss)@f$
     * @return @f$\Gamma(H\to ss)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htotautau 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to tau leptons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to tau leptons 
 * @f$\Gamma(H\to \tau^+\tau^-)@f$ in GeV
 *
 */
class Htotautau: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htotautau(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to tau leptons @f$\Gamma(H\to \tau^+\tau^-)@f$
     * @return @f$\Gamma(H\to \tau^+\tau^-)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htomumu 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to muons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to muons
 * @f$\Gamma(H\to \mu^+\mu^-)@f$ in GeV
 *
 */
class Htomumu: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htomumu(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to muons @f$\Gamma(H\to \mu^+\mu^- )@f$
     * @return @f$\Gamma(H\to \mu^+\mu^-)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class HtoWW 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to @f$WW^*@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to @f$WW^*@f$ 
 * @f$\Gamma(H\to WW^*)@f$ in GeV
 *
 */
class HtoWW: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    HtoWW(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to @f$WW^*@f$ @f$\Gamma(H\to WW^*)@f$
     * @return@f$\Gamma(H\to WW^*)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class HtoZZ 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to @f$ZZ^*@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to @f$ZZ^*@f$ 
 * @f$\Gamma(H\to ZZ^*)@f$ in GeV
 *
 */
class HtoZZ: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    HtoZZ(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The the Higgs decay to @f$ZZ^*@f$ @f$\Gamma(H\to ZZ^*)@f$
     * @return @f$\Gamma(H\to ZZ^*)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htogaga 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to two photons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to two photons 
 * @f$\Gamma(H\to \gamma\gamma)@f$ in GeV
 *
 */
class Htogaga: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htogaga(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to two photons @f$\Gamma(H\to \gamma\gamma)@f$
     * @return @f$\Gamma(H\to \gamma\gamma)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class HtoZga
 * @ingroup EW
 * @brief An observable class for the Higgs decay to Z and photon
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to Z and photon
 * @f$\Gamma(H\to Z\gamma)@f$ in GeV
 *
 */
class HtoZga: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    HtoZga(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to Z and photon @f$\Gamma(H\to Z\gamma)@f$
     * @return @f$\Gamma(H\to Z\gamma)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Htogg 
 * @ingroup EW
 * @brief An observable class for the Higgs decay to gg
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay to gg
 * @f$\Gamma(H\to gg)@f$ in GeV
 *
 */
class Htogg: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Htogg(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay to gg @f$\Gamma(H\to gg)@f$
     * @return @f$\Gamma(H\to gg)@f$ in GeV
     */
    double computeThValue();

    
private:


};

/**
 * @class Hwidth 
 * @ingroup EW
 * @brief An observable class for the Higgs decay width
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs decay width
 * @f$\Gamma(H)@f$ in GeV
 *
 */
class Hwidth: public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Hwidth(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief The Higgs decay width
     * @return @f$\Gamma(H@f$ in GeV
     */
    double computeThValue();

    
private:


};



/** 
 * @}
 */

#endif	/* H_DECAYS_H */
