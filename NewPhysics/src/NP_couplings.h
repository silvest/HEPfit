/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPCOUPLINGS_H
#define	NPCOUPLINGS_H

#include "gslpp.h"
#include "StandardModel.h"
#include "NPbase.h"
#include <ThObservable.h>
#include <string.h>
#include <stdexcept>

//-----  Zff couplings observables  ----------

/**
 * @class deltagZveveL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
 *
 */
class deltagZveveL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZveveL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZveveL class.
     */
    virtual ~deltagZveveL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZvmuvmuL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
 *
 */
class deltagZvmuvmuL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZvmuvmuL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZvmuvmuL class.
     */
    virtual ~deltagZvmuvmuL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZvtavtaL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
 *
 */
class deltagZvtavtaL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZvtavtaL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZvtavtaL class.
     */
    virtual ~deltagZvtavtaL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZeeL
 * @brief An observable class for the deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}/g_{SM}@f$.
 *
 */
class deltagZeeL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZeeL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZeeL class.
     */
    virtual ~deltagZeeL();

    /**
     * @brief The deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling @f$\delta g_{Zee}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZeeR
 * @brief An observable class for the deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}/g_{SM}@f$.
 *
 */
class deltagZeeR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZeeR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZeeR class.
     */
    virtual ~deltagZeeR();

    /**
     * @brief The deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling @f$\delta g_{Zee}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZmumuL
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
 *
 */
class deltagZmumuL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZmumuL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZmumuL class.
     */
    virtual ~deltagZmumuL();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZmumuR
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
 *
 */
class deltagZmumuR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZmumuR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZmumuR class.
     */
    virtual ~deltagZmumuR();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZtataL
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
 *
 */
class deltagZtataL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZtataL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZtataL class.
     */
    virtual ~deltagZtataL();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZtataR
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
 *
 */
class deltagZtataR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZtataR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZtataR class.
     */
    virtual ~deltagZtataR();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuL
 * @brief An observable class for the deviation from the SM of the @f$Z u_{L} u_{L}@f$ coupling
 * @f$\delta g_{Zuu}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u_{L} u_{L}@f$ coupling
 * @f$\delta g_{Zuu}^{L}/g_{SM}@f$.
 *
 */
class deltagZuuL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZuuL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZuuL class.
     */
    virtual ~deltagZuuL();

    /**
     * @brief The deviation from the SM of the @f$Z u_{L} u_{L}@f$ coupling @f$\delta g_{Zuu}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuR
 * @brief An observable class for the deviation from the SM of the @f$Z u_{R} u_{R}@f$ coupling
 * @f$\delta g_{Zuu}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u_{R} u_{R}@f$ coupling
 * @f$\delta g_{Zuu}^{R}/g_{SM}@f$.
 *
 */
class deltagZuuR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZuuR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZuuR class.
     */
    virtual ~deltagZuuR();

    /**
     * @brief The deviation from the SM of the @f$Z u_{R} u_{R}@f$ coupling @f$\delta g_{Zuu}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuV
 * @brief An observable class for the deviation from the SM of the @f$Z u u@f$ vector coupling
 * @f$\delta g_{Zuu}^{V}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u u@f$ vector coupling
 * @f$\delta g_{Zuu}^{V}/g_{SM}@f$.
 *
 */
class deltagZuuV : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZuuV(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZuuV class.
     */
    virtual ~deltagZuuV();

    /**
     * @brief The deviation from the SM of the @f$Z uu@f$ vector coupling @f$\delta g_{Zuu}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuA
 * @brief An observable class for the deviation from the SM of the @f$Z uu@f$ axial coupling
 * @f$\delta g_{Zuu}^{A}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z uu@f$ axial coupling
 * @f$\delta g_{Zuu}^{A}/g_{SM}@f$.
 *
 */
class deltagZuuA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZuuA(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZuuA class.
     */
    virtual ~deltagZuuA();

    /**
     * @brief The deviation from the SM of the @f$Z uu@f$ axial coupling @f$\delta g_{Zuu}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZccL
 * @brief An observable class for the deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
 *
 */
class deltagZccL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZccL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZccL class.
     */
    virtual ~deltagZccL();

    /**
     * @brief The deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZccR
 * @brief An observable class for the deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
 *
 */
class deltagZccR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZccR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZccR class.
     */
    virtual ~deltagZccR();

    /**
     * @brief The deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttL
 * @brief An observable class for the deviation from the SM of the @f$Z t_{L} t_{L}@f$ coupling
 * @f$\delta g_{Ztt}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t_{L} t_{L}@f$ coupling
 * @f$\delta g_{Ztt}^{L}/g_{SM}@f$.
 *
 */
class deltagZttL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZttL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZttL class.
     */
    virtual ~deltagZttL();

    /**
     * @brief The deviation from the SM of the @f$Z t_{L} t_{L}@f$ coupling @f$\delta g_{Ztt}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttR
 * @brief An observable class for the deviation from the SM of the @f$Z t_{R} t_{R}@f$ coupling
 * @f$\delta g_{Ztt}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t_{R} t_{R}@f$ coupling
 * @f$\delta g_{Ztt}^{R}/g_{SM}@f$.
 *
 */
class deltagZttR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZttR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZttR class.
     */
    virtual ~deltagZttR();

    /**
     * @brief The deviation from the SM of the @f$Z t_{R} t_{R}@f$ coupling @f$\delta g_{Ztt}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttV
 * @brief An observable class for the deviation from the SM of the @f$Z t t@f$ vector coupling
 * @f$\delta g_{Ztt}^{V}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t t@f$ vector coupling
 * @f$\delta g_{Ztt}^{V}/g_{SM}@f$.
 *
 */
class deltagZttV : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZttV(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZttV class.
     */
    virtual ~deltagZttV();

    /**
     * @brief The deviation from the SM of the @f$Z t t@f$ vector coupling @f$\delta g_{Ztt}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttA
 * @brief An observable class for the deviation from the SM of the @f$Z t t@f$ axial coupling
 * @f$\delta g_{Ztt}^{A}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t t@f$ axial coupling
 * @f$\delta g_{Ztt}^{A}/g_{SM}@f$.
 *
 */
class deltagZttA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZttA(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZttA class.
     */
    virtual ~deltagZttA();

    /**
     * @brief The deviation from the SM of the @f$Z t t@f$ axial coupling @f$\delta g_{Ztt}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZddL
 * @brief An observable class for the deviation from the SM of the @f$Z d_{L} d_{L}@f$ coupling
 * @f$\delta g_{Zdd}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z d_{L} d_{L}@f$ coupling
 * @f$\delta g_{Zdd}^{L}/g_{SM}@f$.
 *
 */
class deltagZddL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZddL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZddL class.
     */
    virtual ~deltagZddL();

    /**
     * @brief The deviation from the SM of the @f$Z d_{L} d_{L}@f$ coupling @f$\delta g_{Zdd}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddR
 * @brief An observable class for the deviation from the SM of the @f$Z d_{R} d_{R}@f$ coupling
 * @f$\delta g_{Zdd}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z d_{R} d_{R}@f$ coupling
 * @f$\delta g_{Zdd}^{R}/g_{SM}@f$.
 *
 */
class deltagZddR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZddR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZddR class.
     */
    virtual ~deltagZddR();

    /**
     * @brief The deviation from the SM of the @f$Z d_{R} d_{R}@f$ coupling @f$\delta g_{Zdd}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddV
 * @brief An observable class for the deviation from the SM of the @f$Z dd@f$ vector coupling
 * @f$\delta g_{Zdd}^{V}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z dd@f$ vector coupling
 * @f$\delta g_{Zdd}^{V}/g_{SM}@f$.
 *
 */
class deltagZddV : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZddV(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZddV class.
     */
    virtual ~deltagZddV();

    /**
     * @brief The deviation from the SM of the @f$Z dd@f$ vector coupling @f$\delta g_{Zdd}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddA
 * @brief An observable class for the deviation from the SM of the @f$Z dd@f$ axial coupling
 * @f$\delta g_{Zdd}^{A}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z dd@f$ axial coupling
 * @f$\delta g_{Zdd}^{A}/g_{SM}@f$.
 *
 */
class deltagZddA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZddA(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZuuA class.
     */
    virtual ~deltagZddA();

    /**
     * @brief The deviation from the SM of the @f$Z dd@f$ axial coupling @f$\delta g_{Zdd}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZssL
 * @brief An observable class for the deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}/g_{SM}@f$.
 *
 */
class deltagZssL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZssL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZssL class.
     */
    virtual ~deltagZssL();

    /**
     * @brief The deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling @f$\delta g_{Zss}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZssR
 * @brief An observable class for the deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}/g_{SM}@f$.
 *
 */
class deltagZssR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZssR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZssR class.
     */
    virtual ~deltagZssR();

    /**
     * @brief The deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling @f$\delta g_{Zss}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZbbL
 * @brief An observable class for the deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
 *
 */
class deltagZbbL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZbbL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZbbL class.
     */
    virtual ~deltagZbbL();

    /**
     * @brief The deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZbbR
 * @brief An observable class for the deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
 *
 */
class deltagZbbR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagZbbR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagZbbR class.
     */
    virtual ~deltagZbbR();

    /**
     * @brief The deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  Wff couplings observables  ----------

/**
 * @class deltaUWeve
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta U_{We\nu}^{L}/U_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta U_{We\nu}^{L}/U_{SM}@f$.
 *
 */
class deltaUWeve : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaUWeve(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaUWeve class.
     */
    virtual ~deltaUWeve();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta U_{We\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{We\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaUWmuvmu
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$.
 *
 */
class deltaUWmuvmu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaUWmuvmu(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaUWmuvmu class.
     */
    virtual ~deltaUWmuvmu();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaUWtavta
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$.
 *
 */
class deltaUWtavta : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaUWtavta(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaUWtavta class.
     */
    virtual ~deltaUWtavta();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVudL
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling
 * @f$\delta V_{Wud}^{L}/V_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling
 * @f$\delta V_{Wud}^{L}/V_{SM}@f$.
 *
 */
class deltaVudL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVudL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVudL class.
     */
    virtual ~deltaVudL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling @f$\delta V_{Wud}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wud}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVudR
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{u}_{R} d_{R}@f$ coupling
 * @f$\delta V_{Wud}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{u}_{R} d_{R}@f$ coupling
 * @f$\delta V_{Wud}^{R}@f$.
 *
 */
class deltaVudR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVudR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVudR class.
     */
    virtual ~deltaVudR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{u}_{R} d_{R}@f$ coupling @f$\delta V_{Wud}^{R}@f$.
     * @return @f$\delta V_{Wud}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVcsL
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling
 * @f$\delta V_{Wcs}^{L}/V_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling
 * @f$\delta V_{Wcs}^{L}/V_{SM}@f$.
 *
 */
class deltaVcsL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVcsL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVcsL class.
     */
    virtual ~deltaVcsL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling @f$\delta V_{Wcs}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wcs}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVcsR
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{c}_{R} s_{R}@f$ coupling
 * @f$\delta V_{Wcs}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{c}_{R} s_{R}@f$ coupling
 * @f$\delta V_{Wcs}^{R}@f$.
 *
 */
class deltaVcsR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVcsR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVcsR class.
     */
    virtual ~deltaVcsR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{c}_{R} s_{R}@f$ coupling @f$\delta V_{Wcs}^{R}@f$.
     * @return @f$\delta V_{Wcs}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVtbL
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling
 * @f$\delta V_{Wtb}^{L}/V_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling
 * @f$\delta V_{Wtb}^{L}/V_{SM}@f$.
 *
 */
class deltaVtbL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVtbL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVtbL class.
     */
    virtual ~deltaVtbL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling @f$\delta V_{Wtb}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wtb}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVtbR
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{t}_{R} b_{R}@f$ coupling
 * @f$\delta V_{Wtb}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{t}_{R} b_{R}@f$ coupling
 * @f$\delta V_{Wtb}^{R}@f$.
 *
 */
class deltaVtbR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaVtbR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaVtbR class.
     */
    virtual ~deltaVtbR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{t}_{R} b_{R}@f$ coupling @f$\delta V_{Wtb}^{R}@f$.
     * @return @f$\delta V_{Wtb}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


//-----  Hff couplings observables  ----------

/**
 * @class deltagHee
 * @brief An observable class for the deviation from the SM of the @f$H ee@f$ coupling
 * @f$\delta g_{Hee}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H ee@f$ coupling
 * @f$\delta g_{Hee}/g_{SM}@f$.
 *
 */
class deltagHee : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHee(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHee class.
     */
    virtual ~deltagHee();

    /**
     * @brief The deviation from the SM of the @f$H ee@f$ coupling @f$\delta g_{Hee}/g_{SM}@f$.
     * @return @f$\delta g_{Hee}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHmumu
 * @brief An observable class for the deviation from the SM of the @f$H \mu \mu@f$ coupling
 * @f$\delta g_{H\mu\mu}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H \mu \mu@f$ coupling
 * @f$\delta g_{H\mu\mu}/g_{SM}@f$.
 *
 */
class deltagHmumu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHmumu(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHmumu();

    /**
     * @brief The deviation from the SM of the @f$H \mu \mu@f$ coupling @f$\delta g_{H\mu\mu}/g_{SM}@f$.
     * @return @f$\delta g_{H\mu\mu}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHmumueff
 * @brief An observable class for the effective @f$H \mu\mu@f$ coupling 
 * @f$g_{H\mu\mu}^{Eff}@f$, defined from the square root of @f$\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H \mu\mu@f$ coupling
 * @f$g_{H\mu\mu}^{Eff}@f$, defined from the square root of @f$\Gamma_{H\mu\mu}/\Gamma_{H\mu\mu}^{SM}@f$.
 *
 */
class gHmumueff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHmumueff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHmumueff class.
     */
    virtual ~gHmumueff();

    /**
     * @brief The effective @f$H \mu\mu@f$ coupling
     * @return @f$g_{H\mu\mu}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHtata
 * @brief An observable class for the deviation from the SM of the @f$H \tau \tau@f$ coupling
 * @f$\delta g_{H\tau\tau}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H \tau \tau@f$ coupling
 * @f$\delta g_{H\tau\tau}/g_{SM}@f$.
 *
 */
class deltagHtata : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHtata(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHtata();

    /**
     * @brief The deviation from the SM of the @f$H \tau \tau@f$ coupling @f$\delta g_{H\tau\tau}/g_{SM}@f$.
     * @return @f$\delta g_{H\tau\tau}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHtataeff
 * @brief An observable class for the effective @f$H \tau \tau@f$ coupling 
 * @f$g_{H\mu\mu}^{Eff}@f$, defined from the square root of @f$\Gamma_{H\tau \tau}/\Gamma_{H\tau \tau}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H \tau \tau@f$ coupling
 * @f$g_{H\tau \tau}^{Eff}@f$, defined from the square root of @f$\Gamma_{H\tau \tau}/\Gamma_{H\tau \tau}^{SM}@f$.
 *
 */
class gHtataeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHtataeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHtataeff class.
     */
    virtual ~gHtataeff();

    /**
     * @brief The effective @f$H \tau \tau@f$ coupling
     * @return @f$g_{H\tau \tau}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHuu
 * @brief An observable class for the deviation from the SM of the @f$H uu@f$ coupling
 * @f$\delta g_{Huu}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H uu@f$ coupling
 * @f$\delta g_{Huu}/g_{SM}@f$.
 *
 */
class deltagHuu : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHuu(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHuu();

    /**
     * @brief The deviation from the SM of the @f$H uu@f$ coupling @f$\delta g_{Huu}/g_{SM}@f$.
     * @return @f$\delta g_{Huu}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHcc
 * @brief An observable class for the deviation from the SM of the @f$H cc@f$ coupling
 * @f$\delta g_{Hcc}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H cc@f$ coupling
 * @f$\delta g_{Hcc}/g_{SM}@f$.
 *
 */
class deltagHcc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHcc(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHcc();

    /**
     * @brief The deviation from the SM of the @f$H cc@f$ coupling @f$\delta g_{Hcc}/g_{SM}@f$.
     * @return @f$\delta g_{Hcc}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHcceff
 * @brief An observable class for the effective @f$H cc@f$ coupling 
 * @f$g_{Hcc}^{Eff}@f$, defined from the square root of @f$\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H cc@f$ coupling
 * @f$g_{Hcc}^{Eff}@f$, defined from the square root of @f$\Gamma_{Hcc}/\Gamma_{Hcc}^{SM}@f$.
 *
 */
class gHcceff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHcceff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHcceff class.
     */
    virtual ~gHcceff();

    /**
     * @brief The effective @f$H cc@f$ coupling
     * @return @f$g_{Hcc}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHtt
 * @brief An observable class for the deviation from the SM of the @f$H tt@f$ coupling
 * @f$\delta g_{Htt}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H tt@f$ coupling
 * @f$\delta g_{Htt}/g_{SM}@f$.
 *
 */
class deltagHtt : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHtt(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHtt();

    /**
     * @brief The deviation from the SM of the @f$H tt@f$ coupling @f$\delta g_{Htt}/g_{SM}@f$.
     * @return @f$\delta g_{Htt}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHdd
 * @brief An observable class for the deviation from the SM of the @f$H dd@f$ coupling
 * @f$\delta g_{Hdd}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H dd@f$ coupling
 * @f$\delta g_{Hdd}/g_{SM}@f$.
 *
 */
class deltagHdd : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHdd(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHdd class.
     */
    virtual ~deltagHdd();

    /**
     * @brief The deviation from the SM of the @f$H dd@f$ coupling @f$\delta g_{Hdd}/g_{SM}@f$.
     * @return @f$\delta g_{Hdd}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHss
 * @brief An observable class for the deviation from the SM of the @f$H ss@f$ coupling
 * @f$\delta g_{Hss}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H ss@f$ coupling
 * @f$\delta g_{Hss}/g_{SM}@f$.
 *
 */
class deltagHss : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHss(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHss class.
     */
    virtual ~deltagHss();

    /**
     * @brief The deviation from the SM of the @f$H ss@f$ coupling @f$\delta g_{Hss}/g_{SM}@f$.
     * @return @f$\delta g_{Hss}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHbb
 * @brief An observable class for the deviation from the SM of the @f$H b b@f$ coupling
 * @f$\delta g_{Hbb}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H bb@f$ coupling
 * @f$\delta g_{Hbb}/g_{SM}@f$.
 *
 */
class deltagHbb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHbb(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHbb();

    /**
     * @brief The deviation from the SM of the @f$H bb@f$ coupling @f$\delta g_{Hbb}/g_{SM}@f$.
     * @return @f$\delta g_{Hbb}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHbbeff
 * @brief An observable class for the effective @f$H bb@f$ coupling 
 * @f$g_{Hbb}^{Eff}@f$, defined from the square root of @f$\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H bb@f$ coupling
 * @f$g_{Hbb}^{Eff}@f$, defined from the square root of @f$\Gamma_{Hbb}/\Gamma_{Hbb}^{SM}@f$.
 *
 */
class gHbbeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHbbeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHbbeff class.
     */
    virtual ~gHbbeff();

    /**
     * @brief The effective @f$H bb@f$ coupling
     * @return @f$g_{Hbb}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HGG couplings observables  ----------

/**
 * @class deltagHGG
 * @brief An observable class for the deviation from the SM of the effective @f$H g g@f$ coupling
 * @f$\delta g_{HGG}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H g g@f$ coupling
 * @f$\delta g_{HGG}/g_{SM}@f$.
 *
 */
class deltagHGG : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHGG(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHGG class.
     */
    virtual ~deltagHGG();

    /**
     * @brief The deviation from the SM of the effective @f$H g g@f$ coupling @f$\delta g_{HGG}/g_{SM}@f$.
     * @return @f$\delta g_{HGG}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHGGeff
 * @brief An observable class for the effective @f$H GG@f$ coupling 
 * @f$g_{HGG}^{Eff}@f$, defined from the square root of @f$\Gamma_{HGG}/\Gamma_{HGG}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H GG@f$ coupling
 * @f$g_{HGG}^{Eff}@f$, defined from the square root of @f$\Gamma_{HGG}/\Gamma_{HGG}^{SM}@f$.
 *
 */
class gHGGeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHGGeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHGGeff class.
     */
    virtual ~gHGGeff();

    /**
     * @brief The effective @f$H GG@f$ coupling
     * @return @f$g_{HGG}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HZZ couplings observables  ----------

/**
 * @class deltagHZZ
 * @brief An observable class for the deviation from the SM of the @f$H Z Z@f$ coupling
 * @f$\delta g_{HZZ}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H Z Z@f$ coupling
 * @f$\delta g_{HZZ}/g_{SM}@f$.
 *
 */
class deltagHZZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHZZ(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHZZ class.
     */
    virtual ~deltagHZZ();

    /**
     * @brief The deviation from the SM of the effective @f$H Z Z@f$ coupling @f$\delta g_{HZZ}/g_{SM}@f$.
     * @return @f$\delta g_{HZZ}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHZZeff
 * @brief An observable class for the effective @f$H ZZ@f$ coupling 
 * @f$g_{HZZ}^{Eff}@f$, defined from the square root of @f$\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H ZZ@f$ coupling
 * @f$g_{HZZ}^{Eff}@f$, defined from the square root of @f$\Gamma_{HZZ}/\Gamma_{HZZ}^{SM}@f$.
 *
 */
class gHZZeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZZeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZZeff class.
     */
    virtual ~gHZZeff();

    /**
     * @brief The effective @f$H ZZ@f$ coupling
     * @return @f$g_{HZZ}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHZZ1
 * @brief An observable class for the non-SM coupling @f$H Z_{\mu\nu} Z^{\mu\nu}@f$
 * @f$g_{HZZ}^{(1)}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the non-SM coupling @f$H Z_{\mu\nu} Z^{\mu\nu}@f$
 * @f$g_{HZZ}^{(1)}@f$.
 *
 */
class gHZZ1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZZ1(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZZ1 class.
     */
    virtual ~gHZZ1();

    /**
     * @brief The non-SM coupling @f$H Z_{\mu\nu} Z^{\mu\nu}@f$ @f$g_{HZZ}^{(1)}@f$.
     * @return @f$g_{HZZ}^{(1)}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHZZ2
 * @brief An observable class for the non-SM coupling @f$H Z_{\nu} \partial_\mu Z^{\mu\nu}@f$
 * @f$g_{HZZ}^{(2)}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the non-SM coupling @f$H Z_{\nu} \partial_\mu Z^{\mu\nu}@f$
 * @f$g_{HZZ}^{(2)}@f$.
 *
 */
class gHZZ2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZZ2(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZZ2 class.
     */
    virtual ~gHZZ2();

    /**
     * @brief The non-SM coupling @f$H Z_{\nu} \partial_\mu Z^{\mu\nu}@f$ @f$g_{HZZ}^{(2)}@f$.
     * @return @f$g_{HZZ}^{(2)}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HAA couplings observables  ----------

/**
 * @class deltagHAA
 * @brief An observable class for the deviation from the SM of the effective @f$H \gamma \gamma@f$ coupling
 * @f$\delta g_{HAA}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H \gamma \gamma@f$ coupling
 * @f$\delta g_{HAA}/g_{SM}@f$.
 *
 */
class deltagHAA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHAA(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHAA class.
     */
    virtual ~deltagHAA();

    /**
     * @brief The deviation from the SM of the effective @f$H \gamma \gamma@f$ coupling @f$\delta g_{HAA}/g_{SM}@f$.
     * @return @f$\delta g_{HAA}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class gHAAeff
 * @brief An observable class for the effective @f$H AA@f$ coupling 
 * @f$g_{HAA}^{Eff}@f$, defined from the square root of @f$\Gamma_{HAA}/\Gamma_{HAA}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H AA@f$ coupling
 * @f$g_{HAA}^{Eff}@f$, defined from the square root of @f$\Gamma_{HAA}/\Gamma_{HAA}^{SM}@f$.
 *
 */
class gHAAeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHAAeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHAAeff class.
     */
    virtual ~gHAAeff();

    /**
     * @brief The effective @f$H AA@f$ coupling
     * @return @f$g_{HAA}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HZA couplings observables  ----------

/**
 * @class deltagHZA
 * @brief An observable class for the deviation from the SM of the effective @f$H Z \gamma@f$ coupling
 * @f$\delta g_{HZA}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H Z \gamma@f$ coupling
 * @f$\delta g_{HZA}/g_{SM}@f$.
 *
 */
class deltagHZA : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHZA(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHZA class.
     */
    virtual ~deltagHZA();

    /**
     * @brief The deviation from the SM of the effective @f$H Z \gamma@f$ coupling @f$\delta g_{HZA}/g_{SM}@f$.
     * @return @f$\delta g_{HZA}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHZAeff
 * @brief An observable class for the effective @f$H ZA@f$ coupling 
 * @f$g_{HZA}^{Eff}@f$, defined from the square root of @f$\Gamma_{HZA}/\Gamma_{HZA}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H ZA@f$ coupling
 * @f$g_{HZA}^{Eff}@f$, defined from the square root of @f$\Gamma_{HZA}/\Gamma_{HZA}^{SM}@f$.
 *
 */
class gHZAeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZAeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZAeff class.
     */
    virtual ~gHZAeff();

    /**
     * @brief The effective @f$H ZA@f$ coupling
     * @return @f$g_{HZA}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHZA2
 * @brief An observable class for the non-SM coupling @f$H Z_\nu \partial_\mu F^{\mu\nu}@f$
 * @f$g_{HZA}^{(2)}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the non-SM coupling @f$H Z_\nu \partial_\mu F^{\mu\nu}@f$
 * @f$g_{HZA}^{(2)}@f$.
 *
 */
class gHZA2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZA2(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZA2 class.
     */
    virtual ~gHZA2();

    /**
     * @brief The non-SM coupling @f$H Z_\nu \partial_\mu F^{\mu\nu}@f$ @f$g_{HZA}^{(2)}@f$.
     * @return @f$g_{HZA}^{(2)}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HWW couplings observables  ----------

/**
 * @class deltagHWW
 * @brief An observable class for the deviation from the SM of the @f$H W W@f$ coupling
 * @f$\delta g_{HWW}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H W W@f$ coupling
 * @f$\delta g_{HWW}/g_{SM}@f$.
 *
 */
class deltagHWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagHWW(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagHWW class.
     */
    virtual ~deltagHWW();

    /**
     * @brief The deviation from the SM of the effective @f$H W W@f$ coupling @f$\delta g_{HWW}/g_{SM}@f$.
     * @return @f$\delta g_{HWW}/g_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHWWeff
 * @brief An observable class for the effective @f$H WW@f$ coupling 
 * @f$g_{HWW}^{Eff}@f$, defined from the square root of @f$\Gamma_{HWW}/\Gamma_{HWW}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H WW@f$ coupling
 * @f$g_{HWW}^{Eff}@f$, defined from the square root of @f$\Gamma_{HWW}/\Gamma_{HWW}^{SM}@f$.
 *
 */
class gHWWeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWWeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHWWeff class.
     */
    virtual ~gHWWeff();

    /**
     * @brief The effective @f$H WW@f$ coupling
     * @return @f$g_{HWW}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHWW1
 * @brief An observable class for the non-SM coupling @f$H W_{\mu\nu} W^{\mu\nu}@f$
 * @f$g_{HWW}^{(1)}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the non-SM coupling @f$H W_{\mu\nu} W^{\mu\nu}@f$
 * @f$g_{HWW}^{(1)}@f$.
 *
 */
class gHWW1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWW1(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHWW1 class.
     */
    virtual ~gHWW1();

    /**
     * @brief The non-SM coupling @f$H W_{\mu\nu} W^{\mu\nu}@f$ @f$g_{HWW}^{(1)}@f$.
     * @return @f$g_{HWW}^{(1)}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHWW2
 * @brief An observable class for the non-SM coupling @f$H W_{\nu}^+ \partial_\mu W^{-\mu\nu}@f$
 * @f$g_{HWW}^{(2)}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the non-SM coupling @f$H W_{\nu}^+ \partial_\mu W^{-\mu\nu}@f$
 * @f$g_{HWW}^{(2)}@f$.
 *
 */
class gHWW2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWW2(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHWW2 class.
     */
    virtual ~gHWW2();

    /**
     * @brief The non-SM coupling @f$H W_{\nu}^+ \partial_\mu W^{-\mu\nu}@f$ @f$g_{HWW}^{(2)}@f$.
     * @return @f$g_{HWW}^{(2)}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  Other couplings observables  ----------

/**
 * @class gHWZeff
 * @brief An observable class for the ratio of the effective @f$H WW@f$ and @f$H ZZ@f$ couplings 
 * @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the effective @f$H WW@f$ and @f$H ZZ@f$ couplings 
 * @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$.
 *
 */
class gHWZeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWZeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHWZeff class.
     */
    virtual ~gHWZeff();

    /**
     * @brief The ratio of the effective @f$H WW@f$ and @f$H ZZ@f$ couplings
     * @return @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class gHbWeff
 * @brief An observable class for the ratio of the effective @f$H bb@f$ and @f$H WW@f$ couplings 
 * @f$g_{Hbb}^{Eff}/g_{HWW}^{Eff}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the effective @f$H bb@f$ and @f$H WW@f$ couplings 
 * @f$g_{Hbb}^{Eff}/g_{HWW}^{Eff}@f$.
 *
 */
class gHbWeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHbWeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHbWeff class.
     */
    virtual ~gHbWeff();

    /**
     * @brief The ratio of the effective @f$H bb@f$ and @f$H WW@f$ couplings
     * @return @f$g_{Hbb}^{Eff}/g_{HWW}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class gHtaWeff
 * @brief An observable class for the ratio of the effective @f$H \tau\tau@f$ and @f$H WW@f$ couplings 
 * @f$g_{H\tau\tau}^{Eff}/g_{HWW}^{Eff}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the effective @f$H \tau\tau@f$ and @f$H WW@f$ couplings 
 * @f$g_{H\tau\tau}^{Eff}/g_{HWW}^{Eff}@f$.
 *
 */
class gHtaWeff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHtaWeff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHtaWeff class.
     */
    virtual ~gHtaWeff();

    /**
     * @brief The ratio of the effective @f$H \tau\tau@f$ and @f$H WW@f$ couplings
     * @return @f$g_{H\tau\tau}^{Eff}/g_{HWW}^{Eff}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HHH couplings observables  ----------

/**
 * @class deltalHHH
 * @brief An observable class for the deviation from the SM of the @f$H H H@f$ coupling
 * @f$\delta \lambda_{H^3}/\lambda_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H H H@f$ coupling
 * @f$\delta \lambda_{H^3}/\lambda_{SM}@f$.
 *
 */
class deltalHHH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltalHHH(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltalHHH class.
     */
    virtual ~deltalHHH();

    /**
     * @brief The deviation from the SM of the effective @f$H H H@f$ coupling @f$\delta \lambda_{H^3}/\lambda_{SM}@f$.
     * @return @f$\delta \lambda_{H^3}/\lambda_{SM}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  VVV couplings observables  ----------

// See aTGC in EW




//-----  Basic interactions of the so-called Higgs basis  ----------

/**
 * @class deltaytHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta y_t@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta y_t@f$.
 *
 */
class deltaytHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaytHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaytHB class.
     */
    virtual ~deltaytHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * @return @f$\delta y_t@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaybHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta y_b@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta y_b@f$.
 *
 */
class deltaybHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaybHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaybHB class.
     */
    virtual ~deltaybHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * @return @f$\delta y_b@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaytauHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta y_\tau@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta y_\tau@f$.
 *
 */
class deltaytauHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaytauHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaytauHB class.
     */
    virtual ~deltaytauHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * @return @f$\delta y_\tau@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaycHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta y_c@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta y_c@f$.
 *
 */
class deltaycHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaycHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaycHB class.
     */
    virtual ~deltaycHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * @return @f$\delta y_c@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltaymuHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta y_\mu@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta y_\mu@f$.
 *
 */
class deltaymuHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaymuHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaymuHB class.
     */
    virtual ~deltaymuHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * @return @f$\delta y_\mu@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltacZHB
 * @brief An observable class for the Higgs-basis coupling @f$\delta c_z@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\delta c_z@f$.
 *
 */
class deltacZHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltacZHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltacZHB class.
     */
    virtual ~deltacZHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * @return @f$\delta c_z@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class cZBoxHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{z\Box}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{z\Box}@f$.
 *
 */
class cZBoxHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cZBoxHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the cZBoxHB class.
     */
    virtual ~cZBoxHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * @return @f$c_{z\Box}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class cZZHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{zz}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{zz}@f$.
 *
 */
class cZZHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cZZHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the cZZHB class.
     */
    virtual ~cZZHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * @return @f$c_{zz}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class cZgaHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{z\gamma}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{z\gamma}@f$.
 *
 */
class cZgaHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cZgaHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the cZgaHB class.
     */
    virtual ~cZgaHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * @return @f$c_{z\gamma}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class cgagaHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{\gamma\gamma}@f$.
 *
 */
class cgagaHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cgagaHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the cgagaHB class.
     */
    virtual ~cgagaHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * @return @f$c_{\gamma\gamma}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class cggHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{gg}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{gg}@f$.
 *
 */
class cggHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cggHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the cggHB class.
     */
    virtual ~cggHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * @return @f$c_{gg}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class lambzHB
 * @brief An observable class for the Higgs-basis coupling @f$\lambda_{z}@f$.
 * (See LHCHXSWG-INT-2015-001 document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$\lambda_{z}@f$.
 *
 */
class lambzHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    lambzHB(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the lambzHB class.
     */
    virtual ~lambzHB();

    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * @return @f$\lambda_{z}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


//-----  Other useful observables to work with new physics  ----------

/**
 * @class deltaMW
 * @brief An observable class for the deviation from the SM of the @f$W@f$ mass
 * @f$\delta M_{W}/M_{W}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM of the @f$W@f$ mass
 * @f$\delta M_{W}/M_{W}@f$.
 *
 */
class deltaMW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaMW(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaMW class.
     */
    virtual ~deltaMW();

    /**
     * @brief The deviation from the SM of the @f$W@f$ mass.
     * @return @f$\delta M_{W}/M_{W}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


#endif	/* NPCOUPLINGS_H */

