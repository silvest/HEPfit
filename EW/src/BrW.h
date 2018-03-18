/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BRW_H
#define	BRW_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class BrWlepton 
 * @ingroup EW 
 * @brief An observable class for the branching ratio of the @f$W@f$ boson into 
 * a given lepton + neutrino flavour.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the branching ratio of the @f$W@f$ boson 
 * into a given lepton + neutrino flavour, @f$BR(W\to \ell \nu)@f$. 
 * Assumes lepton-flavour universality.
 */
class BrWlepton : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    BrWlepton(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The the branching ratio of the @f$W@f$ boson
     * into a given lepton + neutrino flavour, @f$BR(W\to \ell \nu)@f$.
     * @return @f$BR(W\to \ell \nu)@f$
     */
    double computeThValue();

};

/**
 * @class BrWelectron 
 * @ingroup EW 
 * @brief An observable class for the branching ratio of the @f$W@f$ boson into 
 * an electron + neutrino.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the branching ratio of the @f$W@f$ boson 
 * into an electron + neutrino, @f$BR(W\to e \nu)@f$. 
 * Lepton-flavour universality not assumed.
 */
class BrWelectron : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    BrWelectron(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The the branching ratio of the @f$W@f$ boson
     * into an electron + neutrino, @f$BR(W\to e \nu)@f$.
     * @return @f$BR(W\to e \nu)@f$
     */
    double computeThValue();

};

/**
 * @class BrWmuon 
 * @ingroup EW 
 * @brief An observable class for the branching ratio of the @f$W@f$ boson into 
 * a muon + neutrino.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the branching ratio of the @f$W@f$ boson 
 * into a muon + neutrino, @f$BR(W\to \mu \nu)@f$. 
 * Lepton-flavour universality not assumed.
 */
class BrWmuon : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    BrWmuon(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The the branching ratio of the @f$W@f$ boson
     * into a muon + neutrino, @f$BR(W\to \mu \nu)@f$.
     * @return @f$BR(W\to \mu \nu)@f$
     */
    double computeThValue();

};

/**
 * @class BrWtau 
 * @ingroup EW 
 * @brief An observable class for the branching ratio of the @f$W@f$ boson into 
 * a tau + neutrino.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the branching ratio of the @f$W@f$ boson 
 * into a tau + neutrino, @f$BR(W\to \tau \nu)@f$. 
 * Lepton-flavour universality not assumed.
 */
class BrWtau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    BrWtau(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The the branching ratio of the @f$W@f$ boson
     * into a tau + neutrino, @f$BR(W\to \tau \nu)@f$.
     * @return @f$BR(W\to \tau \nu)@f$
     */
    double computeThValue();

};

/**
 * @class BrWhadrons 
 * @ingroup EW 
 * @brief An observable class for the branching ratio of the @f$W@f$ boson into 
 * hadrons.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the branching ratio of the @f$W@f$ boson 
 * into hadrons, @f$BR(W\to had)@f$. 
 * Lepton-flavour universality not assumed.
 */
class BrWhadrons : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    BrWhadrons(const StandardModel& SM_i) 
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The the branching ratio of the @f$W@f$ boson
     * into hadrons, @f$BR(W\to hadrons)@f$.
     * @return @f$BR(W\to had)@f$
     */
    double computeThValue();

};

#endif	/* BRW_H */

