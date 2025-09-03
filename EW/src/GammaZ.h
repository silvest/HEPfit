/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAZ_H
#define	GAMMAZ_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class GammaZ
 * @ingroup EW 
 * @brief An observable class for the total decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the total decay width of the @f$Z@f$
 * boson,
 * @f[
 * \Gamma_Z
 * = 3\,\Gamma_\nu + \Gamma_{e} + \Gamma_{\mu} + \Gamma_{\tau} + \Gamma_h\,,
 * @f]
 * where @f$\Gamma_h=\sum_{q\neq t}\Gamma_q@f$ is the total hadronic width.
 *
 * @sa EW_NPZff::GammaZ() and the detailed description of EW class
 * for the inclusion of new physics contribution
 * 
 */
class GammaZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZ(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$,
     * in units of GeV.
     * @return @f$\Gamma_Z@f$ in units of GeV
     */
    double computeThValue();

    
private:


};




/**
 * @class GammaZee
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZee : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZee(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to ee)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to ee)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZmumu
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZmumu(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to \mu\mu)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to \mu\mu)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZtautau
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZtautau : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZtautau(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to \tau\tau)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to \tau\tau)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZuu
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZuu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZuu(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to uu)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to uu)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZcc
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZcc(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to cc)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to cc)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZdd
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZdd : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZdd(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to dd)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to dd)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZss
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZss : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZss(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to ss)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to ss)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZbb
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZbb(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to bb)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to bb)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};


/**
 * @class GammaZinv
 * @ingroup EW 
 * @brief An observable class for the partial decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the partial decay width of the @f$Z@f$
 * boson 
 * 
 */
class GammaZinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZinv(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The partial decay width of the @f$Z@f$ boson, @f$\Gamma(Z\to \nu\nu)@f$,
     * in units of GeV.
     * @return @f$\Gamma(Z\to \nu\nu)@f$ in units of GeV
     */
    double computeThValue();

    
private:


};




#endif	/* GAMMAZ_H */

