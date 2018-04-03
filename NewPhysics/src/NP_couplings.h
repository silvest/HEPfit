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
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}@f$.
     * @return @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZvmuvmuL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}@f$.
     * @return @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZvtavtaL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}@f$.
     * @return @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZeeL
 * @brief An observable class for the deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling @f$\delta g_{Zee}^{L}@f$.
     * @return @f$\delta g_{Zee}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZeeR
 * @brief An observable class for the deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling @f$\delta g_{Zee}^{R}@f$.
     * @return @f$\delta g_{Zee}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZmumuL
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling @f$\delta g_{Z\mu\mu}^{L}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZmumuR
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling @f$\delta g_{Z\mu\mu}^{R}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZtataL
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling @f$\delta g_{Z\tau\tau}^{L}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZtataR
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling @f$\delta g_{Z\tau\tau}^{R}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuL
 * @brief An observable class for the deviation from the SM of the @f$Z u_{L} u_{L}@f$ coupling
 * @f$\delta g_{Zuu}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u_{L} u_{L}@f$ coupling
 * @f$\delta g_{Zuu}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z u_{L} u_{L}@f$ coupling @f$\delta g_{Zuu}^{L}@f$.
     * @return @f$\delta g_{Zuu}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuR
 * @brief An observable class for the deviation from the SM of the @f$Z u_{R} u_{R}@f$ coupling
 * @f$\delta g_{Zuu}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u_{R} u_{R}@f$ coupling
 * @f$\delta g_{Zuu}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z u_{R} u_{R}@f$ coupling @f$\delta g_{Zuu}^{R}@f$.
     * @return @f$\delta g_{Zuu}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuV
 * @brief An observable class for the deviation from the SM of the @f$Z u u@f$ vector coupling
 * @f$\delta g_{Zuu}^{V}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z u u@f$ vector coupling
 * @f$\delta g_{Zuu}^{V}@f$.
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
     * @brief The deviation from the SM of the @f$Z uu@f$ vector coupling @f$\delta g_{Zuu}^{V}@f$.
     * @return @f$\delta g_{Zuu}^{V}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZuuA
 * @brief An observable class for the deviation from the SM of the @f$Z uu@f$ axial coupling
 * @f$\delta g_{Zuu}^{A}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z uu@f$ axial coupling
 * @f$\delta g_{Zuu}^{A}@f$.
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
     * @brief The deviation from the SM of the @f$Z uu@f$ axial coupling @f$\delta g_{Zuu}^{A}@f$.
     * @return @f$\delta g_{Zuu}^{A}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZccL
 * @brief An observable class for the deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling @f$\delta g_{Zcc}^{L}@f$.
     * @return @f$\delta g_{Zcc}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZccR
 * @brief An observable class for the deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling @f$\delta g_{Zcc}^{R}@f$.
     * @return @f$\delta g_{Zcc}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttL
 * @brief An observable class for the deviation from the SM of the @f$Z t_{L} t_{L}@f$ coupling
 * @f$\delta g_{Ztt}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t_{L} t_{L}@f$ coupling
 * @f$\delta g_{Ztt}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z t_{L} t_{L}@f$ coupling @f$\delta g_{Ztt}^{L}@f$.
     * @return @f$\delta g_{Ztt}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttR
 * @brief An observable class for the deviation from the SM of the @f$Z t_{R} t_{R}@f$ coupling
 * @f$\delta g_{Ztt}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t_{R} t_{R}@f$ coupling
 * @f$\delta g_{Ztt}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z t_{R} t_{R}@f$ coupling @f$\delta g_{Ztt}^{R}@f$.
     * @return @f$\delta g_{Ztt}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttV
 * @brief An observable class for the deviation from the SM of the @f$Z t t@f$ vector coupling
 * @f$\delta g_{Ztt}^{V}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t t@f$ vector coupling
 * @f$\delta g_{Ztt}^{V}@f$.
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
     * @brief The deviation from the SM of the @f$Z t t@f$ vector coupling @f$\delta g_{Ztt}^{V}@f$.
     * @return @f$\delta g_{Ztt}^{V}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZttA
 * @brief An observable class for the deviation from the SM of the @f$Z t t@f$ axial coupling
 * @f$\delta g_{Ztt}^{A}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z t t@f$ axial coupling
 * @f$\delta g_{Ztt}^{A}@f$.
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
     * @brief The deviation from the SM of the @f$Z t t@f$ axial coupling @f$\delta g_{Ztt}^{A}@f$.
     * @return @f$\delta g_{Ztt}^{A}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagZddL
 * @brief An observable class for the deviation from the SM of the @f$Z d_{L} d_{L}@f$ coupling
 * @f$\delta g_{Zdd}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z d_{L} d_{L}@f$ coupling
 * @f$\delta g_{Zdd}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z d_{L} d_{L}@f$ coupling @f$\delta g_{Zdd}^{L}@f$.
     * @return @f$\delta g_{Zdd}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddR
 * @brief An observable class for the deviation from the SM of the @f$Z d_{R} d_{R}@f$ coupling
 * @f$\delta g_{Zdd}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z d_{R} d_{R}@f$ coupling
 * @f$\delta g_{Zdd}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z d_{R} d_{R}@f$ coupling @f$\delta g_{Zdd}^{R}@f$.
     * @return @f$\delta g_{Zdd}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddV
 * @brief An observable class for the deviation from the SM of the @f$Z dd@f$ vector coupling
 * @f$\delta g_{Zdd}^{V}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z dd@f$ vector coupling
 * @f$\delta g_{Zdd}^{V}@f$.
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
     * @brief The deviation from the SM of the @f$Z dd@f$ vector coupling @f$\delta g_{Zdd}^{V}@f$.
     * @return @f$\delta g_{Zdd}^{V}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZddA
 * @brief An observable class for the deviation from the SM of the @f$Z dd@f$ axial coupling
 * @f$\delta g_{Zdd}^{A}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z dd@f$ axial coupling
 * @f$\delta g_{Zdd}^{A}@f$.
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
     * @brief The deviation from the SM of the @f$Z dd@f$ axial coupling @f$\delta g_{Zdd}^{A}@f$.
     * @return @f$\delta g_{Zdd}^{A}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZssL
 * @brief An observable class for the deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling @f$\delta g_{Zss}^{L}@f$.
     * @return @f$\delta g_{Zss}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZssR
 * @brief An observable class for the deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling @f$\delta g_{Zss}^{R}@f$.
     * @return @f$\delta g_{Zss}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZbbL
 * @brief An observable class for the deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling @f$\delta g_{Zbb}^{L}@f$.
     * @return @f$\delta g_{Zbb}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagZbbR
 * @brief An observable class for the deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}@f$.
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
     * @brief The deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling @f$\delta g_{Zbb}^{R}@f$.
     * @return @f$\delta g_{Zbb}^{R}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  Wff couplings observables  ----------

/**
 * @class deltaUWeve
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta U_{We\nu}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta U_{We\nu}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta U_{We\nu}^{L}@f$.
     * @return @f$\delta U_{We\nu}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaUWmuvmu
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta U_{W\mu\nu}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta U_{W\mu\nu}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta U_{W\mu\nu}^{L}@f$.
     * @return @f$\delta U_{W\mu\nu}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaUWtavta
 * @brief An observable class for the deviation from the SM of the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta U_{W\tau\nu}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta U_{W\tau\nu}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta U_{W\tau\nu}^{L}@f$.
     * @return @f$\delta U_{W\tau\nu}^{L}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltaVudL
 * @brief An observable class for the deviation from the SM of the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling
 * @f$\delta V_{Wud}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling
 * @f$\delta V_{Wud}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling @f$\delta V_{Wud}^{L}@f$.
     * @return @f$\delta V_{Wud}^{L}@f$
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
 * @f$\delta V_{Wcs}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling
 * @f$\delta V_{Wcs}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling @f$\delta V_{Wcs}^{L}@f$.
     * @return @f$\delta V_{Wcs}^{L}@f$
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
 * @f$\delta V_{Wtb}^{L}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling
 * @f$\delta V_{Wtb}^{L}@f$.
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
     * @brief The deviation from the SM of the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling @f$\delta V_{Wtb}^{L}@f$.
     * @return @f$\delta V_{Wtb}^{L}@f$
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
 * @f$\delta g_{Hee}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H ee@f$ coupling
 * @f$\delta g_{Hee}@f$.
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
     * @brief The deviation from the SM of the @f$H ee@f$ coupling @f$\delta g_{Hee}@f$.
     * @return @f$\delta g_{Hee}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHmumu
 * @brief An observable class for the deviation from the SM of the @f$H \mu \mu@f$ coupling
 * @f$\delta g_{H\mu\mu}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H \mu \mu@f$ coupling
 * @f$\delta g_{H\mu\mu}@f$.
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
     * @brief The deviation from the SM of the @f$H \mu \mu@f$ coupling @f$\delta g_{H\mu\mu}@f$.
     * @return @f$\delta g_{H\mu\mu}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHtata
 * @brief An observable class for the deviation from the SM of the @f$H \tau \tau@f$ coupling
 * @f$\delta g_{H\tau\tau}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H \tau \tau@f$ coupling
 * @f$\delta g_{H\tau\tau}@f$.
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
     * @brief The deviation from the SM of the @f$H \tau \tau@f$ coupling @f$\delta g_{H\tau\tau}@f$.
     * @return @f$\delta g_{H\tau\tau}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHuu
 * @brief An observable class for the deviation from the SM of the @f$H uu@f$ coupling
 * @f$\delta g_{Huu}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H uu@f$ coupling
 * @f$\delta g_{Huu}@f$.
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
     * @brief The deviation from the SM of the @f$H uu@f$ coupling @f$\delta g_{Huu}@f$.
     * @return @f$\delta g_{Huu}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHcc
 * @brief An observable class for the deviation from the SM of the @f$H cc@f$ coupling
 * @f$\delta g_{Hcc}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H cc@f$ coupling
 * @f$\delta g_{Hcc}@f$.
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
     * @brief The deviation from the SM of the @f$H cc@f$ coupling @f$\delta g_{Hcc}@f$.
     * @return @f$\delta g_{Hcc}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHtt
 * @brief An observable class for the deviation from the SM of the @f$H tt@f$ coupling
 * @f$\delta g_{Htt}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H tt@f$ coupling
 * @f$\delta g_{Htt}@f$.
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
     * @brief The deviation from the SM of the @f$H tt@f$ coupling @f$\delta g_{Htt}@f$.
     * @return @f$\delta g_{Htt}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHdd
 * @brief An observable class for the deviation from the SM of the @f$H dd@f$ coupling
 * @f$\delta g_{Hdd}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H dd@f$ coupling
 * @f$\delta g_{Hdd}@f$.
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
     * @brief The deviation from the SM of the @f$H dd@f$ coupling @f$\delta g_{Hdd}@f$.
     * @return @f$\delta g_{Hdd}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class deltagHss
 * @brief An observable class for the deviation from the SM of the @f$H ss@f$ coupling
 * @f$\delta g_{Hss}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H ss@f$ coupling
 * @f$\delta g_{Hss}@f$.
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
     * @brief The deviation from the SM of the @f$H ss@f$ coupling @f$\delta g_{Hss}@f$.
     * @return @f$\delta g_{Hss}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class deltagHbb
 * @brief An observable class for the deviation from the SM of the @f$H b b@f$ coupling
 * @f$\delta g_{Hbb}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$H bb@f$ coupling
 * @f$\delta g_{Hbb}@f$.
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
     * @brief The deviation from the SM of the @f$H bb@f$ coupling @f$\delta g_{Hbb}@f$.
     * @return @f$\delta g_{Hbb}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HGG couplings observables  ----------

/**
 * @class deltagHGG
 * @brief An observable class for the deviation from the SM of the effective @f$H g g@f$ coupling
 * @f$\delta g_{HGG}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H g g@f$ coupling
 * @f$\delta g_{HGG}@f$.
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
     * @brief The deviation from the SM of the effective @f$H g g@f$ coupling @f$\delta g_{HGG}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HZZ couplings observables  ----------

/**
 * @class deltagHZZ
 * @brief An observable class for the deviation from the SM of the @f$H Z Z@f$ coupling
 * @f$\delta g_{HZZ}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H Z Z@f$ coupling
 * @f$\delta g_{HZZ}@f$.
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
     * @brief The deviation from the SM of the effective @f$H Z Z@f$ coupling @f$\delta g_{HZZ}@f$.
     * @return @f$\delta g_{HZZ}@f$
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
 * @f$\delta g_{HAA}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H \gamma \gamma@f$ coupling
 * @f$\delta g_{HAA}@f$.
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
     * @brief The deviation from the SM of the effective @f$H \gamma \gamma@f$ coupling @f$\delta g_{HAA}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  HZA couplings observables  ----------

/**
 * @class deltagHZA
 * @brief An observable class for the deviation from the SM of the effective @f$H Z \gamma@f$ coupling
 * @f$\delta g_{HZA}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H Z \gamma@f$ coupling
 * @f$\delta g_{HZA}@f$.
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
     * @brief The deviation from the SM of the effective @f$H Z \gamma@f$ coupling @f$\delta g_{HZA}@f$.
     * @return @f$\delta g_{HZA}@f$
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
 * @f$\delta g_{HWW}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H W W@f$ coupling
 * @f$\delta g_{HWW}@f$.
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
     * @brief The deviation from the SM of the effective @f$H W W@f$ coupling @f$\delta g_{HWW}@f$.
     * @return @f$\delta g_{HWW}@f$
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

//-----  HHH couplings observables  ----------

/**
 * @class deltalHHH
 * @brief An observable class for the deviation from the SM of the @f$H H H@f$ coupling
 * @f$\delta \lambda_{H^3}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the effective @f$H H H@f$ coupling
 * @f$\delta \lambda_{H^3}@f$.
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
     * @brief The deviation from the SM of the effective @f$H H H@f$ coupling @f$\delta \lambda_{H^3}@f$.
     * @return @f$\delta \lambda_{H^3}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

//-----  VVV couplings observables  ----------

// See aTGC in EW

#endif	/* NPCOUPLINGS_H */

