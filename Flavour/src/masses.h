#ifndef MASSES_H
#define MASSES_H

#include "ThObservable.h"
class StandardModel;

/**
 * @class up_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{u}(2\mathrm{GeV})@f$, the running mass of the up quark at 2 GeV
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{u}(2\mathrm{GeV})@f$.
 */
class up_mass : public ThObservable {
public:

    /**
    * @brief Constructor declaration.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    up_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{u}(2\mathrm{GeV})@f$ 
     */
    double computeThValue();

};

/**
 * @class down_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{d}(2\mathrm{GeV})@f$, the running mass of the down quark at 2 GeV
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{d}(2\mathrm{GeV})@f$.
 */
class down_mass : public ThObservable {
public:

    /**
    * @brief Constructor declaration.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    down_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{d}(2\mathrm{GeV})@f$ 
     */
    virtual double computeThValue();

};

/**
 * @class strange_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{s}(2\mathrm{GeV})@f$, the running mass of the strange quark at 2 GeV
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{s}(2\mathrm{GeV})@f$.
 */
class strange_mass : public ThObservable {
public:

    /**
    * @brief Constructor declaration.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    strange_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{s}(2\mathrm{GeV})@f$ 
     */
    virtual double computeThValue();

};

/**
 * @class charm_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{c}(m_{c})@f$, the running mass of the charm quark at the charm mass scale @f$m_{c}@f$ (in GeV) 
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{c}(m_{c})@f$.
 */
class charm_mass : public ThObservable {
public:

    /**
    * @brief Constructor declaration.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    charm_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{c}(m_{c})@f$ 
     */
    virtual double computeThValue();

};

/**
 * @class bottom_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{b}(m_{b})@f$, the running mass of the bottom quark at the bottom mass scale @f$m_{b}@f$ (in GeV) 
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{b}(m_{b})@f$.
 */
class bottom_mass : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    bottom_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{b}(m_{b})@f$ 
     */
    virtual double computeThValue();

};

/**
 * @class top_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{t}(m_{t})@f$, the running mass of the top quark
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{t}(m_{t})@f$.
 */
class top_mass : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    top_mass(const StandardModel& SM_i) ;

    /**
     * 
     * @return value of @f$m_{t}(m_{t})@f$ 
     */
    virtual double computeThValue();
};

/**
 * @class mtpole
 * @ingroup Flavour
 * @brief A class for @f$m_{t}^{\text{pole}}@f$, the pole mass of the top quark
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{t}^{\text{pole}}@f$.
 */
class mtpole : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    mtpole(const StandardModel& SM_i) ;

    /**
     * 
     * @return value of @f$m_{t}^{\text{pole}}@f$ 
     */
    virtual double computeThValue();
};

/**
 * @class electron_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{e}@f$, the mass of the electron 
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{e}@f$.
 */
class electron_mass : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    electron_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{e}(2\mathrm{GeV})@f$ 
     */
    virtual double computeThValue();
};

/**
 * @class muon_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{\mu}@f$, the mass of the muon
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{\mu}@f$.
 */
class muon_mass : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    muon_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{\mu}@f$ 
     */
    virtual double computeThValue();
 };

/**
 * @class tau_mass
 * @ingroup Flavour
 * @brief A class for @f$m_{\tau}@f$, the running mass of the tau lepton
 * @author HEPfit Collaboration
 * @details This class is used to retrieve the value of
 * @f$m_{\tau}@f$.
 */
class tau_mass : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    tau_mass(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$m_{\tau}@f$ 
     */
    virtual double computeThValue();
};

#endif /* MASSES_H */