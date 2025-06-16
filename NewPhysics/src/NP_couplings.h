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

#include <ThObservable.h>
#include <string.h>
#include <stdexcept>

class NPbase;

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
    deltagZveveL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZveveL class.
     */
    virtual ~deltagZveveL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZvmuvmuL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZvmuvmuL class.
     */
    virtual ~deltagZvmuvmuL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$
     */
    double computeThValue();
    
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZvtavtaL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZvtavtaL class.
     */
    virtual ~deltagZvtavtaL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZeeL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZeeL class.
     */
    virtual ~deltagZeeL();

    /**
     * @brief The deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling @f$\delta g_{Zee}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZeeR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZeeR class.
     */
    virtual ~deltagZeeR();

    /**
     * @brief The deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling @f$\delta g_{Zee}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZmumuL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZmumuL class.
     */
    virtual ~deltagZmumuL();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZmumuR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZmumuR class.
     */
    virtual ~deltagZmumuR();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZtataL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZtataL class.
     */
    virtual ~deltagZtataL();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZtataR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZtataR class.
     */
    virtual ~deltagZtataR();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZuuL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZuuL class.
     */
    virtual ~deltagZuuL();

    /**
     * @brief The deviation from the SM of the @f$Z u_{L} u_{L}@f$ coupling @f$\delta g_{Zuu}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZuuR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZuuR class.
     */
    virtual ~deltagZuuR();

    /**
     * @brief The deviation from the SM of the @f$Z u_{R} u_{R}@f$ coupling @f$\delta g_{Zuu}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZuuV(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZuuV class.
     */
    virtual ~deltagZuuV();

    /**
     * @brief The deviation from the SM of the @f$Z uu@f$ vector coupling @f$\delta g_{Zuu}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZuuA(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZuuA class.
     */
    virtual ~deltagZuuA();

    /**
     * @brief The deviation from the SM of the @f$Z uu@f$ axial coupling @f$\delta g_{Zuu}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Zuu}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZccL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZccL class.
     */
    virtual ~deltagZccL();

    /**
     * @brief The deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZccR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZccR class.
     */
    virtual ~deltagZccR();

    /**
     * @brief The deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZttL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZttL class.
     */
    virtual ~deltagZttL();

    /**
     * @brief The deviation from the SM of the @f$Z t_{L} t_{L}@f$ coupling @f$\delta g_{Ztt}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZttR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZttR class.
     */
    virtual ~deltagZttR();

    /**
     * @brief The deviation from the SM of the @f$Z t_{R} t_{R}@f$ coupling @f$\delta g_{Ztt}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZttV(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZttV class.
     */
    virtual ~deltagZttV();

    /**
     * @brief The deviation from the SM of the @f$Z t t@f$ vector coupling @f$\delta g_{Ztt}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZttA(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZttA class.
     */
    virtual ~deltagZttA();

    /**
     * @brief The deviation from the SM of the @f$Z t t@f$ axial coupling @f$\delta g_{Ztt}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Ztt}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZddL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZddL class.
     */
    virtual ~deltagZddL();

    /**
     * @brief The deviation from the SM of the @f$Z d_{L} d_{L}@f$ coupling @f$\delta g_{Zdd}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZddR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZddR class.
     */
    virtual ~deltagZddR();

    /**
     * @brief The deviation from the SM of the @f$Z d_{R} d_{R}@f$ coupling @f$\delta g_{Zdd}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZddV(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZddV class.
     */
    virtual ~deltagZddV();

    /**
     * @brief The deviation from the SM of the @f$Z dd@f$ vector coupling @f$\delta g_{Zdd}^{V}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{V}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZddA(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZuuA class.
     */
    virtual ~deltagZddA();

    /**
     * @brief The deviation from the SM of the @f$Z dd@f$ axial coupling @f$\delta g_{Zdd}^{A}/g_{SM}@f$.
     * @return @f$\delta g_{Zdd}^{A}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZssL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZssL class.
     */
    virtual ~deltagZssL();

    /**
     * @brief The deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling @f$\delta g_{Zss}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZssR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZssR class.
     */
    virtual ~deltagZssR();

    /**
     * @brief The deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling @f$\delta g_{Zss}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZbbL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZbbL class.
     */
    virtual ~deltagZbbL();

    /**
     * @brief The deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagZbbR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagZbbR class.
     */
    virtual ~deltagZbbR();

    /**
     * @brief The deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
};


//-----  Zff EFFECTIVE couplings observables: relative corrections (derived from Af and Gamma(Z->ff)  ----------


/**
 * @class deltagEffZveveL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZveveL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZveveL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZveveL class.
     */
    virtual ~deltagEffZveveL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{e}\nu^{e}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZvmuvmuL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZvmuvmuL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZvmuvmuL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZvmuvmuL class.
     */
    virtual ~deltagEffZvmuvmuL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\mu}\nu^{\mu}}^{L}/g_{SM}@f$
     */
    double computeThValue();
    
private:
    const NPbase * myNPbase;

};


/**
 * @class deltagEffZvtavtaL
 * @brief An observable class for the deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling
 * @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZvtavtaL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZvtavtaL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZvtavtaL class.
     */
    virtual ~deltagEffZvtavtaL();

    /**
     * @brief The deviation from the SM of the @f$Z \nu^{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\nu^{\tau}\nu^{\tau}}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};


/**
 * @class deltagEffZeeL
 * @brief An observable class for the deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{L} e_{L}@f$ coupling
 * @f$\delta g_{Zee}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZeeL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZeeL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZeeL class.
     */
    virtual ~deltagEffZeeL();

    /**
     * @brief The deviation from the SM of the @f$Z e_{L} e_{L}@f$ coupling @f$\delta g_{Zee}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZeeR
 * @brief An observable class for the deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z e_{R} e_{R}@f$ coupling
 * @f$\delta g_{Zee}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZeeR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZeeR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZeeR class.
     */
    virtual ~deltagEffZeeR();

    /**
     * @brief The deviation from the SM of the @f$Z e_{R} e_{R}@f$ coupling @f$\delta g_{Zee}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zee}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZmumuL
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{L} \mu_{L}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZmumuL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZmumuL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZmumuL class.
     */
    virtual ~deltagEffZmumuL();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{L} \mu_{L}@f$ coupling @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZmumuR
 * @brief An observable class for the deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \mu_{R} \mu_{R}@f$ coupling
 * @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZmumuR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZmumuR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZmumuR class.
     */
    virtual ~deltagEffZmumuR();

    /**
     * @brief The deviation from the SM of the @f$Z \mu_{R} \mu_{R}@f$ coupling @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\mu\mu}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZtataL
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{L} \tau_{L}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZtataL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZtataL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZtataL class.
     */
    virtual ~deltagEffZtataL();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{L} \tau_{L}@f$ coupling @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZtataR
 * @brief An observable class for the deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z \tau_{R} \tau_{R}@f$ coupling
 * @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZtataR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZtataR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZtataR class.
     */
    virtual ~deltagEffZtataR();

    /**
     * @brief The deviation from the SM of the @f$Z \tau_{R} \tau_{R}@f$ coupling @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Z\tau\tau}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZccL
 * @brief An observable class for the deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{L} c_{L}@f$ coupling
 * @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZccL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZccL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZccL class.
     */
    virtual ~deltagEffZccL();

    /**
     * @brief The deviation from the SM of the @f$Z c_{L} c_{L}@f$ coupling @f$\delta g_{Zcc}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZccR
 * @brief An observable class for the deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z c_{R} c_{R}@f$ coupling
 * @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZccR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZccR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZccR class.
     */
    virtual ~deltagEffZccR();

    /**
     * @brief The deviation from the SM of the @f$Z c_{R} c_{R}@f$ coupling @f$\delta g_{Zcc}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zcc}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};


/**
 * @class deltagEffZssL
 * @brief An observable class for the deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{L} s_{L}@f$ coupling
 * @f$\delta g_{Zss}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZssL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZssL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZssL class.
     */
    virtual ~deltagEffZssL();

    /**
     * @brief The deviation from the SM of the @f$Z s_{L} s_{L}@f$ coupling @f$\delta g_{Zss}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZssR
 * @brief An observable class for the deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z s_{R} s_{R}@f$ coupling
 * @f$\delta g_{Zss}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZssR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZssR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZssR class.
     */
    virtual ~deltagEffZssR();

    /**
     * @brief The deviation from the SM of the @f$Z s_{R} s_{R}@f$ coupling @f$\delta g_{Zss}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zss}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZbbL
 * @brief An observable class for the deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{L} b_{L}@f$ coupling
 * @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
 *
 */
class deltagEffZbbL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZbbL(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZbbL class.
     */
    virtual ~deltagEffZbbL();

    /**
     * @brief The deviation from the SM of the @f$Z b_{L} b_{L}@f$ coupling @f$\delta g_{Zbb}^{L}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{L}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

};

/**
 * @class deltagEffZbbR
 * @brief An observable class for the deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the deviation from the SM on the @f$Z b_{R} b_{R}@f$ coupling
 * @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
 *
 */
class deltagEffZbbR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltagEffZbbR(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltagEffZbbR class.
     */
    virtual ~deltagEffZbbR();

    /**
     * @brief The deviation from the SM of the @f$Z b_{R} b_{R}@f$ coupling @f$\delta g_{Zbb}^{R}/g_{SM}@f$.
     * @return @f$\delta g_{Zbb}^{R}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;

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
    deltaUWeve(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaUWeve class.
     */
    virtual ~deltaUWeve();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{e}_{L} \nu^{e}_{L}@f$ coupling @f$\delta U_{We\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{We\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaUWmuvmu(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaUWmuvmu class.
     */
    virtual ~deltaUWmuvmu();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{\mu}_{L} \nu^{\mu}_{L}@f$ coupling @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{W\mu\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaUWtavta(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaUWtavta class.
     */
    virtual ~deltaUWtavta();

    /**
     * @brief The deviation from the SM of the @f$W^{-} \bar{\tau}_{L} \nu^{\tau}_{L}@f$ coupling @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$.
     * @return @f$\delta U_{W\tau\nu}^{L}/U_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVudL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVudL class.
     */
    virtual ~deltaVudL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{u}_{L} d_{L}@f$ coupling @f$\delta V_{Wud}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wud}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVudR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVudR class.
     */
    virtual ~deltaVudR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{u}_{R} d_{R}@f$ coupling @f$\delta V_{Wud}^{R}@f$.
     * @return @f$\delta V_{Wud}^{R}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVcsL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVcsL class.
     */
    virtual ~deltaVcsL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{c}_{L} s_{L}@f$ coupling @f$\delta V_{Wcs}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wcs}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVcsR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVcsR class.
     */
    virtual ~deltaVcsR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{c}_{R} s_{R}@f$ coupling @f$\delta V_{Wcs}^{R}@f$.
     * @return @f$\delta V_{Wcs}^{R}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVtbL(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVtbL class.
     */
    virtual ~deltaVtbL();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{t}_{L} b_{L}@f$ coupling @f$\delta V_{Wtb}^{L}/V_{SM}@f$.
     * @return @f$\delta V_{Wtb}^{L}/V_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaVtbR(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaVtbR class.
     */
    virtual ~deltaVtbR();

    /**
     * @brief The deviation from the SM of the @f$W^{+} \bar{t}_{R} b_{R}@f$ coupling @f$\delta V_{Wtb}^{R}@f$.
     * @return @f$\delta V_{Wtb}^{R}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHee(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHee class.
     */
    virtual ~deltagHee();

    /**
     * @brief The deviation from the SM of the @f$H ee@f$ coupling @f$\delta g_{Hee}/g_{SM}@f$.
     * @return @f$\delta g_{Hee}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHmumu(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHmumu();

    /**
     * @brief The deviation from the SM of the @f$H \mu \mu@f$ coupling @f$\delta g_{H\mu\mu}/g_{SM}@f$.
     * @return @f$\delta g_{H\mu\mu}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHtata(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHtata();

    /**
     * @brief The deviation from the SM of the @f$H \tau \tau@f$ coupling @f$\delta g_{H\tau\tau}/g_{SM}@f$.
     * @return @f$\delta g_{H\tau\tau}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHuu(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHuu();

    /**
     * @brief The deviation from the SM of the @f$H uu@f$ coupling @f$\delta g_{Huu}/g_{SM}@f$.
     * @return @f$\delta g_{Huu}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHcc(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHcc();

    /**
     * @brief The deviation from the SM of the @f$H cc@f$ coupling @f$\delta g_{Hcc}/g_{SM}@f$.
     * @return @f$\delta g_{Hcc}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHtt(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHtt();

    /**
     * @brief The deviation from the SM of the @f$H tt@f$ coupling @f$\delta g_{Htt}/g_{SM}@f$.
     * @return @f$\delta g_{Htt}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHdd(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHdd class.
     */
    virtual ~deltagHdd();

    /**
     * @brief The deviation from the SM of the @f$H dd@f$ coupling @f$\delta g_{Hdd}/g_{SM}@f$.
     * @return @f$\delta g_{Hdd}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHss(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHss class.
     */
    virtual ~deltagHss();

    /**
     * @brief The deviation from the SM of the @f$H ss@f$ coupling @f$\delta g_{Hss}/g_{SM}@f$.
     * @return @f$\delta g_{Hss}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHbb(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHtata class.
     */
    virtual ~deltagHbb();

    /**
     * @brief The deviation from the SM of the @f$H bb@f$ coupling @f$\delta g_{Hbb}/g_{SM}@f$.
     * @return @f$\delta g_{Hbb}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHGG(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHGG class.
     */
    virtual ~deltagHGG();

    /**
     * @brief The deviation from the SM of the effective @f$H g g@f$ coupling @f$\delta g_{HGG}/g_{SM}@f$.
     * @return @f$\delta g_{HGG}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHZZ(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHZZ class.
     */
    virtual ~deltagHZZ();

    /**
     * @brief The deviation from the SM of the effective @f$H Z Z@f$ coupling @f$\delta g_{HZZ}/g_{SM}@f$.
     * @return @f$\delta g_{HZZ}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
};

/**
 * @class gHZZ4feff
 * @brief An observable class for the effective @f$H ZZ@f$ coupling, defined from 4f final sates 
 * @f$g_{HZZ,4f}^{Eff}@f$, defined from the square root of @f$\Gamma_{H4f, NC}/\Gamma_{H4f, NC}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H ZZ@f$ coupling
 * @f$g_{HZZ,4f}^{Eff}@f$, defined from the square root of @f$\Gamma_{H4f, NC}/\Gamma_{H4f, NC}^{SM}@f$.
 *
 */
class gHZZ4feff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHZZ4feff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHZZeff class.
     */
    virtual ~gHZZ4feff();

    /**
     * @brief The effective @f$H ZZ@f$ coupling
     * @return @f$g_{HZZ,4f}^{Eff}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
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
    gHZZ1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHZZ1 class.
     */
    virtual ~gHZZ1();

    /**
     * @brief The non-SM coupling @f$H Z_{\mu\nu} Z^{\mu\nu}@f$ @f$g_{HZZ}^{(1)}@f$.
     * @return @f$g_{HZZ}^{(1)}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    gHZZ2(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHZZ2 class.
     */
    virtual ~gHZZ2();

    /**
     * @brief The non-SM coupling @f$H Z_{\nu} \partial_\mu Z^{\mu\nu}@f$ @f$g_{HZZ}^{(2)}@f$.
     * @return @f$g_{HZZ}^{(2)}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHAA(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHAA class.
     */
    virtual ~deltagHAA();

    /**
     * @brief The deviation from the SM of the effective @f$H \gamma \gamma@f$ coupling @f$\delta g_{HAA}/g_{SM}@f$.
     * @return @f$\delta g_{HAA}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    deltagHZA(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHZA class.
     */
    virtual ~deltagHZA();

    /**
     * @brief The deviation from the SM of the effective @f$H Z \gamma@f$ coupling @f$\delta g_{HZA}/g_{SM}@f$.
     * @return @f$\delta g_{HZA}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
    gHZA2(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHZA2 class.
     */
    virtual ~gHZA2();

    /**
     * @brief The non-SM coupling @f$H Z_\nu \partial_\mu F^{\mu\nu}@f$ @f$g_{HZA}^{(2)}@f$.
     * @return @f$g_{HZA}^{(2)}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltagHWW(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltagHWW class.
     */
    virtual ~deltagHWW();

    /**
     * @brief The deviation from the SM of the effective @f$H W W@f$ coupling @f$\delta g_{HWW}/g_{SM}@f$.
     * @return @f$\delta g_{HWW}/g_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
};

/**
 * @class gHWW4feff
 * @brief An observable class for the effective @f$H WW@f$ coupling defined from all 4f final states minus pure NC
 * (Still contains a small amount of NC and NC-CC interference)
 * @f$g_{HWW,4f}^{Eff}@f$, defined from the square root of @f$\Gamma_{H4f, CC}/\Gamma_{H4f, CC}^{SM}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective @f$H WW@f$ coupling
 * @f$g_{HWW,4f}^{Eff}@f$, defined from the square root of @f$\Gamma_{H4f, CC}/\Gamma_{H4f, CC}^{SM}@f$.
 *
 */
class gHWW4feff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWW4feff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the gHWWeff class.
     */
    virtual ~gHWW4feff();

    /**
     * @brief The effective @f$H WW@f$ coupling
     * @return @f$g_{HWW,4f}^{Eff}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
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
    gHWW1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHWW1 class.
     */
    virtual ~gHWW1();

    /**
     * @brief The non-SM coupling @f$H W_{\mu\nu} W^{\mu\nu}@f$ @f$g_{HWW}^{(1)}@f$.
     * @return @f$g_{HWW}^{(1)}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    gHWW2(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHWW2 class.
     */
    virtual ~gHWW2();

    /**
     * @brief The non-SM coupling @f$H W_{\nu}^+ \partial_\mu W^{-\mu\nu}@f$ @f$g_{HWW}^{(2)}@f$.
     * @return @f$g_{HWW}^{(2)}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
};

/**
 * @class gHWZSMLin
 * @brief An observable class for the ratio of the SM-like @f$H WW@f$ and @f$H ZZ@f$ couplings 
 * @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$  (linear in new physics effects).
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the effective @f$H WW@f$ and @f$H ZZ@f$ couplings 
 * @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$.
 *
 */
class gHWZSMLin : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    gHWZSMLin(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the gHWZSMLin class.
     */
    virtual ~gHWZSMLin();

    /**
     * @brief The ratio of the SM-like @f$H WW@f$ and @f$H ZZ@f$ couplings
     * @return @f$g_{HWW}^{Eff}/g_{HZZ}^{Eff}@f$ (linear in new physics effects)
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
      
private:
    const NPbase * myNPbase;
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
      
private:
    const NPbase * myNPbase;
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
    deltalHHH(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltalHHH class.
     */
    virtual ~deltalHHH();

    /**
     * @brief The deviation from the SM of the effective @f$H H H@f$ coupling @f$\delta \lambda_{H^3}/\lambda_{SM}@f$.
     * @return @f$\delta \lambda_{H^3}/\lambda_{SM}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
};

//-----  VVV couplings observables  ----------

// See aTGC in EW. Here we define only the Effective couplings used in arXiv: 1708.09079 [hep-ph]

/**
 * @class deltag1ZEff
 * @brief An observable class for the effective anomalous triple gauge coupling
 * @f$\delta g_{1,Z}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective anomalous triple gauge coupling
 * @f$\delta g_{1,Z}^{Eff}@f$.
 *
 */
class deltag1ZEff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltag1ZEff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltag1ZEff class.
     */
    virtual ~deltag1ZEff();

    /**
     * @brief The anomalous triple gauge coupling @f$\delta g_{1,Z}^{Eff}@f$.
     * @return @f$\delta g_{1,Z}^{Eff}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class deltaKgammaEff
 * @brief An observable class for the effective anomalous triple gauge coupling
 * @f$\delta \kappa_{\gamma}^{Eff}@f$ from arXiv: 1708.09079 [hep-ph].
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the effective anomalous triple gauge coupling
 * @f$\delta \kappa_{\gamma}^{Eff}@f$.
 *
 */
class deltaKgammaEff : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    deltaKgammaEff(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the deltaKgammaEff class.
     */
    virtual ~deltaKgammaEff();

    /**
     * @brief The anomalous triple gauge coupling @f$\delta \kappa_{\gamma}^{Eff}@f$.
     * @return @f$\delta \kappa_{\gamma}^{Eff}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


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
    deltaytHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaytHB class.
     */
    virtual ~deltaytHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_t@f$.
     * @return @f$\delta y_t@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaybHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaybHB class.
     */
    virtual ~deltaybHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_b@f$.
     * @return @f$\delta y_b@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaytauHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaytauHB class.
     */
    virtual ~deltaytauHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\tau@f$.
     * @return @f$\delta y_\tau@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaycHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaycHB class.
     */
    virtual ~deltaycHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_c@f$.
     * @return @f$\delta y_c@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltaymuHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltaymuHB class.
     */
    virtual ~deltaymuHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta y_\mu@f$.
     * @return @f$\delta y_\mu@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    deltacZHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltacZHB class.
     */
    virtual ~deltacZHB();

    /**
     * @brief The Higgs-basis coupling @f$\delta c_z@f$.
     * @return @f$\delta c_z@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    cZBoxHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cZBoxHB class.
     */
    virtual ~cZBoxHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{z\Box}@f$.
     * @return @f$c_{z\Box}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    cZZHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cZZHB class.
     */
    virtual ~cZZHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{zz}@f$.
     * @return @f$c_{zz}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    cZgaHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cZgaHB class.
     */
    virtual ~cZgaHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{z\gamma}@f$.
     * @return @f$c_{z\gamma}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    cgagaHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cgagaHB class.
     */
    virtual ~cgagaHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{\gamma\gamma}@f$.
     * @return @f$c_{\gamma\gamma}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    cggHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cggHB class.
     */
    virtual ~cggHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{gg}@f$.
     * @return @f$c_{gg}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
};


/**
 * @class cggEffHB
 * @brief An observable class for the Higgs-basis coupling @f$c_{gg}^{Eff}@f$.
 * (Similar to cgg_HB but including modifications of SM loops.)
 * (See arXiv: 1505.00046 [hep-ph] document.)
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Higgs-basis coupling
 * @f$c_{gg}^{Eff}@f$.
 *
 */
class cggEffHB : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    cggEffHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the cggEffHB class.
     */
    virtual ~cggEffHB();

    /**
     * @brief The Higgs-basis coupling @f$c_{gg}^{Eff}@f$.
     * @return @f$c_{gg}^{Eff}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
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
    lambzHB(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the lambzHB class.
     */
    virtual ~lambzHB();

    /**
     * @brief The Higgs-basis coupling @f$\lambda_{z}@f$.
     * @return @f$\lambda_{z}@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
    const double mu;
};


//-----  Other useful observables to work with new physics  ----------


//-----  Relative correction to W mass  ----------

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
      
private:
    const NPbase * myNPbase;
};


//-----  Relative correction to EM constant  ----------

/**
 * @class deltae
 * @brief An observable class for the relative deviation from the SM of the EM coupling @f$e@f$ 
 * @f$\delta e@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute deviation from the SM of the EM coupling @f$e@f$
 * @f$\delta e@f$.
 *
 */
class deltae : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] lepton a lepton
     */
    deltae(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the deltae class.
     */
    virtual ~deltae();

    /**
     * @brief The absolute deviation from the SM from the SM of the EM coupling @f$e@f$
     * @f$\delta e@f$.
     * @return @f$\delta e@f$
     */
    double computeThValue();    
    
private:
    const NPbase * myNPbase;
    const double mu;
};


//-----  Absolute correction to some EW couplings (factoring e/sc or e/sqrt(2)s  ----------

/**
 * @class delgZlL
 * @brief An observable class for the absolute deviation from the SM of the @f$Z l_{L} l_{L}@f$ coupling
 * @f$\delta g_{Zll}^{L}@f$, factoring out the @f$e/sc@f$ overall coupling.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute deviation from the SM on the @f$Z l_{L} l_{L}@f$ coupling
 * @f$\delta g_{Zll}^{L}@f$.
 *
 */
class delgZlL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] lepton a lepton
     */
    delgZlL(const StandardModel& SM_i, const StandardModel::lepton lepton, const double mu_i);
      
    /**
     * @brief Destructor of the delgZlL class.
     */
    virtual ~delgZlL();

    /**
     * @brief The absolute deviation from the SM of the @f$Z l_{L} l_{L}@f$ coupling
     * @f$\delta g_{Zll}^{L}@f$, factoring out the @f$e/sc@f$ overall coupling.
     * @return @f$\delta g_{Zll}^{L}@f$
     */
    double computeThValue();    
    
private:
    const NPbase * myNPbase;
    StandardModel::lepton lepton;
    const double mu;

};

/**
 * @class delgZlR
 * @brief An observable class for the absolute deviation from the SM of the @f$Z l_{R} l_{R}@f$ coupling
 * @f$\delta g_{Zll}^{R}@f$, factoring out the @f$e/sc@f$ overall coupling.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute deviation from the SM on the @f$Z l_{R} l_{R}@f$ coupling
 * @f$\delta g_{Zll}^{R}@f$.
 *
 */
class delgZlR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] lepton a lepton
     */
    delgZlR(const StandardModel& SM_i, const StandardModel::lepton lepton, const double mu_i);
      
    /**
     * @brief Destructor of the delgZlR class.
     */
    virtual ~delgZlR();

    /**
     * @brief The absolute deviation from the SM of the @f$Z l_{R} l_{R}@f$ coupling
     * @f$\delta g_{Zll}^{R}@f$, factoring out the @f$e/sc@f$ overall coupling.
     * @return @f$\delta g_{Zll}^{R}@f$
     */
    double computeThValue();
    
private:
    const NPbase * myNPbase;
    StandardModel::lepton lepton;
    const double mu;

};

/**
 * @class delgZqL
 * @brief An observable class for the absolute deviation from the SM of the @f$Z q_{L} q_{L}@f$ coupling
 * @f$\delta g_{Zqq}^{L}@f$, factoring out the @f$e/sc@f$ overall coupling.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute deviation from the SM on the @f$Z q_{L} q_{L}@f$ coupling
 * @f$\delta g_{Zqq}^{L}@f$.
 *
 */
class delgZqL : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark a quark
     */
    delgZqL(const StandardModel& SM_i, const StandardModel::quark quark, const double mu_i);
      
    /**
     * @brief Destructor of the delgZqL class.
     */
    virtual ~delgZqL();

    /**
     * @brief The absolute deviation from the SM of the @f$Z q_{L} q_{L}@f$ coupling
     * @f$\delta g_{Zqq}^{L}@f$, factoring out the @f$e/sc@f$ overall coupling.
     * @return @f$\delta g_{Zqq}^{L}@f$
     */
    double computeThValue();    
    
private:
    const NPbase * myNPbase;
    StandardModel::quark quark;
    const double mu;

};

/**
 * @class delgZqR
 * @brief An observable class for the absolute deviation from the SM of the @f$Z q_{R} q_{R}@f$ coupling
 * @f$\delta g_{Zqq}^{R}@f$, factoring out the @f$e/sc@f$ overall coupling.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute deviation from the SM on the @f$Z q_{R} q_{R}@f$ coupling
 * @f$\delta g_{Zqq}^{R}@f$.
 *
 */
class delgZqR : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark a quark
     */
    delgZqR(const StandardModel& SM_i, const StandardModel::quark quark, const double mu_i);
      
    /**
     * @brief Destructor of the delgZqR class.
     */
    virtual ~delgZqR();

    /**
     * @brief The absolute deviation from the SM of the @f$Z q_{R} q_{R}@f$ coupling
     * @f$\delta g_{Zqq}^{R}@f$, factoring out the @f$e/sc@f$ overall coupling.
     * @return @f$\delta g_{Zqq}^{R}@f$
     */
    double computeThValue();
    
private:
    const NPbase * myNPbase;
    StandardModel::quark quark;
    const double mu;

};


//-----  Oblique parameters  ----------

/**
 * @class oblS
 * @brief An observable class for the oblique parameter @f$S@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter @f$S@f$.
 *
 */
class oblS : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    oblS(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~oblS();

    /**
     * @brief The oblique parameter @f$S@f$.
     * @return @f$S@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};

/**
 * @class oblT
 * @brief An observable class for the oblique parameter @f$T@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter @f$T@f$.
 *
 */
class oblT : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    oblT(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblT class.
     */
    virtual ~oblT();

    /**
     * @brief The oblique parameter @f$T@f$.
     * @return @f$T@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class oblW
 * @brief An observable class for the oblique parameter @f$W@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter @f$W@f$.
 *
 */
class oblW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    oblW(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~oblW();

    /**
     * @brief The oblique parameter @f$W@f$.
     * @return @f$W@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class oblY
 * @brief An observable class for the oblique parameter @f$Y@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the oblique parameter @f$Y@f$.
 *
 */
class oblY : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    oblY(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblY class.
     */
    virtual ~oblY();

    /**
     * @brief The oblique parameter @f$Y@f$.
     * @return @f$Y@f$
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


//-----  Combinations of Warsaw basis coefficients constrained by EWPO  ----------

/**
 * @class CEWHL111
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(1)})_{11}@f$.
 *
 */
class CEWHL111 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL111(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL111 class.
     */
    virtual ~CEWHL111();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHL122
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(1)})_{22}@f$.
 *
 */
class CEWHL122 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL122(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL122 class.
     */
    virtual ~CEWHL122();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHL133
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(1)})_{33}@f$.
 *
 */
class CEWHL133 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL133(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL133 class.
     */
    virtual ~CEWHL133();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(1)})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHL311
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(3)})_{11}@f$.
 *
 */
class CEWHL311 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL311(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL311 class.
     */
    virtual ~CEWHL311();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHL322
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(3)})_{22}@f$.
 *
 */
class CEWHL322 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL322(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL322 class.
     */
    virtual ~CEWHL322();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHL333
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HL}^{(3)})_{33}@f$.
 *
 */
class CEWHL333 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHL333(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHL333 class.
     */
    virtual ~CEWHL333();

    /**
     * @brief The combination @f$(\hat{C}_{HL}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HL}^{(3)})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHQ111
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(1)})_{11}@f$.
 *
 */
class CEWHQ111 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ111(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ111 class.
     */
    virtual ~CEWHQ111();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(1)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHQ122
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(1)})_{22}@f$.
 *
 */
class CEWHQ122 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ122(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ122 class.
     */
    virtual ~CEWHQ122();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(1)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHQ133
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(1)})_{33}@f$.
 *
 */
class CEWHQ133 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ133(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ133 class.
     */
    virtual ~CEWHQ133();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(1)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(1)})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHQ311
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(3)})_{11}@f$.
 *
 */
class CEWHQ311 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ311(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ311 class.
     */
    virtual ~CEWHQ311();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(3)})_{11}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHQ322
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(3)})_{22}@f$.
 *
 */
class CEWHQ322 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ322(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ322 class.
     */
    virtual ~CEWHQ322();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(3)})_{22}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};

/**
 * @class CEWHQ333
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(3)})_{33}@f$.
 *
 */
class CEWHQ333 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQ333(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ333 class.
     */
    virtual ~CEWHQ333();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(3)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(3)})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHQd33
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{HQ}^{(d)})_{33}@f$.
 *
 */
class CEWHQd33 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHQd33(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHQ333 class.
     */
    virtual ~CEWHQd33();

    /**
     * @brief The combination @f$(\hat{C}_{HQ}^{(d)})_{33}@f$.
     * @return @f$(\hat{C}_{HQ}^{(d)})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHe11
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{He})_{11}@f$.
 *
 */
class CEWHe11 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHe11(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHe11 class.
     */
    virtual ~CEWHe11();

    /**
     * @brief The combination @f$(\hat{C}_{He})_{11}@f$.
     * @return @f$(\hat{C}_{He})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHe22
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{He})_{22}@f$.
 *
 */
class CEWHe22 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHe22(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHe22 class.
     */
    virtual ~CEWHe22();

    /**
     * @brief The combination @f$(\hat{C}_{He})_{22}@f$.
     * @return @f$(\hat{C}_{He})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHe33
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{He})_{33}@f$.
 *
 */
class CEWHe33 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHe33(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHe33 class.
     */
    virtual ~CEWHe33();

    /**
     * @brief The combination @f$(\hat{C}_{He})_{33}@f$.
     * @return @f$(\hat{C}_{He})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHu11
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hu})_{11}@f$.
 *
 */
class CEWHu11 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHu11(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHu11 class.
     */
    virtual ~CEWHu11();

    /**
     * @brief The combination @f$(\hat{C}_{Hu})_{11}@f$.
     * @return @f$(\hat{C}_{Hu})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHu22
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hu})_{22}@f$.
 *
 */
class CEWHu22 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHu22(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHu22 class.
     */
    virtual ~CEWHu22();

    /**
     * @brief The combination @f$(\hat{C}_{Hu})_{22}@f$.
     * @return @f$(\hat{C}_{Hu})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHu33
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hu})_{33}@f$.
 *
 */
class CEWHu33 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHu33(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHu33 class.
     */
    virtual ~CEWHu33();

    /**
     * @brief The combination @f$(\hat{C}_{Hu})_{33}@f$.
     * @return @f$(\hat{C}_{Hu})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHd11
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hd})_{11}@f$.
 *
 */
class CEWHd11 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHd11(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHd11 class.
     */
    virtual ~CEWHd11();

    /**
     * @brief The combination @f$(\hat{C}_{Hd})_{11}@f$.
     * @return @f$(\hat{C}_{Hd})_{11}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHd22
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hd})_{22}@f$.
 *
 */
class CEWHd22 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHd22(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHd22 class.
     */
    virtual ~CEWHd22();

    /**
     * @brief The combination @f$(\hat{C}_{Hd})_{22}@f$.
     * @return @f$(\hat{C}_{Hd})_{22}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


/**
 * @class CEWHd33
 * @brief An observable class for the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the combinations of coefficients of the Warsaw basis constrained by EWPO
 * @f$(\hat{C}_{Hd})_{33}@f$.
 *
 */
class CEWHd33 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    CEWHd33(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the CEWHd33 class.
     */
    virtual ~CEWHd33();

    /**
     * @brief The combination @f$(\hat{C}_{Hd})_{33}@f$.
     * @return @f$(\hat{C}_{Hd})_{33}@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:


};


//-----  Auxiliary observables  ----------

/**
 * @class AuxObsNP1
 * @brief An observable class for the auxiliary observable AuxObsNP1
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP1.
 *
 */
class AuxObsNP1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP1(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP1();

    /**
     * @brief The auxiliary observable AuxObsNP1.
     * @return AuxObsNP1
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP2
 * @brief An observable class for the auxiliary observable AuxObsNP2
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP2.
 *
 */
class AuxObsNP2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP2(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP2();

    /**
     * @brief The auxiliary observable AuxObsNP2.
     * @return AuxObsNP2
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP3
 * @brief An observable class for the auxiliary observable AuxObsNP3
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP3.
 *
 */
class AuxObsNP3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP3(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP3();

    /**
     * @brief The auxiliary observable AuxObsNP3.
     * @return AuxObsNP3
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP4
 * @brief An observable class for the auxiliary observable AuxObsNP4
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP4.
 *
 */
class AuxObsNP4 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP4(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP4();

    /**
     * @brief The auxiliary observable AuxObsNP4.
     * @return AuxObsNP4
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP5
 * @brief An observable class for the auxiliary observable AuxObsNP5
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP5.
 *
 */
class AuxObsNP5 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP5(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP5();

    /**
     * @brief The auxiliary observable AuxObsNP5.
     * @return AuxObsNP5
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP6
 * @brief An observable class for the auxiliary observable AuxObsNP6
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP6.
 *
 */
class AuxObsNP6 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP6(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP6();

    /**
     * @brief The auxiliary observable AuxObsNP6.
     * @return AuxObsNP6
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP7
 * @brief An observable class for the auxiliary observable AuxObsNP7
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP7.
 *
 */
class AuxObsNP7 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP7(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP7();

    /**
     * @brief The auxiliary observable AuxObsNP7.
     * @return AuxObsNP7
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP8
 * @brief An observable class for the auxiliary observable AuxObsNP8
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP8.
 *
 */
class AuxObsNP8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP8(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP8();

    /**
     * @brief The auxiliary observable AuxObsNP8.
     * @return AuxObsNP8
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP9
 * @brief An observable class for the auxiliary observable AuxObsNP9
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP9.
 *
 */
class AuxObsNP9 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP9(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP9();

    /**
     * @brief The auxiliary observable AuxObsNP9.
     * @return AuxObsNP9
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP10
 * @brief An observable class for the auxiliary observable AuxObsNP10
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP10.
 *
 */
class AuxObsNP10 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP10(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP10();

    /**
     * @brief The auxiliary observable AuxObsNP10.
     * @return AuxObsNP10
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP11
 * @brief An observable class for the auxiliary observable AuxObsNP11
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP11.
 *
 */
class AuxObsNP11 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP11(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP11();

    /**
     * @brief The auxiliary observable AuxObsNP11.
     * @return AuxObsNP11
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP12
 * @brief An observable class for the auxiliary observable AuxObsNP12
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP12.
 *
 */
class AuxObsNP12 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP12(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP12();

    /**
     * @brief The auxiliary observable AuxObsNP12.
     * @return AuxObsNP12
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};



/**
 * @class AuxObsNP13
 * @brief An observable class for the auxiliary observable AuxObsNP13
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP13.
 *
 */
class AuxObsNP13 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP13(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP13();

    /**
     * @brief The auxiliary observable AuxObsNP13.
     * @return AuxObsNP13
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP14
 * @brief An observable class for the auxiliary observable AuxObsNP14
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP14.
 *
 */
class AuxObsNP14 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP14(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP14();

    /**
     * @brief The auxiliary observable AuxObsNP14.
     * @return AuxObsNP14
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP15
 * @brief An observable class for the auxiliary observable AuxObsNP15
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP15.
 *
 */
class AuxObsNP15 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP15(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP15();

    /**
     * @brief The auxiliary observable AuxObsNP15.
     * @return AuxObsNP15
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP16
 * @brief An observable class for the auxiliary observable AuxObsNP16
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP16.
 *
 */
class AuxObsNP16 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP16(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP16();

    /**
     * @brief The auxiliary observable AuxObsNP16.
     * @return AuxObsNP16
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP17
 * @brief An observable class for the auxiliary observable AuxObsNP17
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP17.
 *
 */
class AuxObsNP17 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP17(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP17();

    /**
     * @brief The auxiliary observable AuxObsNP17.
     * @return AuxObsNP17
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP18
 * @brief An observable class for the auxiliary observable AuxObsNP18
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP18.
 *
 */
class AuxObsNP18 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP18(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP18();

    /**
     * @brief The auxiliary observable AuxObsNP18.
     * @return AuxObsNP18
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP19
 * @brief An observable class for the auxiliary observable AuxObsNP19
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP19.
 *
 */
class AuxObsNP19 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP19(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP19();

    /**
     * @brief The auxiliary observable AuxObsNP19.
     * @return AuxObsNP19
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP20
 * @brief An observable class for the auxiliary observable AuxObsNP20
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP20.
 *
 */
class AuxObsNP20 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP20(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP20();

    /**
     * @brief The auxiliary observable AuxObsNP20.
     * @return AuxObsNP20
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP21
 * @brief An observable class for the auxiliary observable AuxObsNP21
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP21.
 *
 */
class AuxObsNP21 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP21(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP21();

    /**
     * @brief The auxiliary observable AuxObsNP21.
     * @return AuxObsNP21
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP22
 * @brief An observable class for the auxiliary observable AuxObsNP22
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP22.
 *
 */
class AuxObsNP22 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP22(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP22();

    /**
     * @brief The auxiliary observable AuxObsNP22.
     * @return AuxObsNP22
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};



/**
 * @class AuxObsNP23
 * @brief An observable class for the auxiliary observable AuxObsNP23
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP23.
 *
 */
class AuxObsNP23 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP23(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP23();

    /**
     * @brief The auxiliary observable AuxObsNP23.
     * @return AuxObsNP23
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP24
 * @brief An observable class for the auxiliary observable AuxObsNP24
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP24.
 *
 */
class AuxObsNP24 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP24(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP24();

    /**
     * @brief The auxiliary observable AuxObsNP24.
     * @return AuxObsNP24
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP25
 * @brief An observable class for the auxiliary observable AuxObsNP25
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP25.
 *
 */
class AuxObsNP25 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP25(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP25();

    /**
     * @brief The auxiliary observable AuxObsNP25.
     * @return AuxObsNP25
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP26
 * @brief An observable class for the auxiliary observable AuxObsNP26
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP26.
 *
 */
class AuxObsNP26 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP26(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP26();

    /**
     * @brief The auxiliary observable AuxObsNP26.
     * @return AuxObsNP26
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP27
 * @brief An observable class for the auxiliary observable AuxObsNP27
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP27.
 *
 */
class AuxObsNP27 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP27(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP27();

    /**
     * @brief The auxiliary observable AuxObsNP27.
     * @return AuxObsNP27
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP28
 * @brief An observable class for the auxiliary observable AuxObsNP28
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP28.
 *
 */
class AuxObsNP28 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP28(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP28();

    /**
     * @brief The auxiliary observable AuxObsNP28.
     * @return AuxObsNP28
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP29
 * @brief An observable class for the auxiliary observable AuxObsNP29
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP29.
 *
 */
class AuxObsNP29 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP29(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP29();

    /**
     * @brief The auxiliary observable AuxObsNP29.
     * @return AuxObsNP29
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};


/**
 * @class AuxObsNP30
 * @brief An observable class for the auxiliary observable AuxObsNP20
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the auxiliary observable AuxObsNP20.
 *
 */
class AuxObsNP30 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AuxObsNP30(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the oblW class.
     */
    virtual ~AuxObsNP30();

    /**
     * @brief The auxiliary observable AuxObsNP30.
     * @return AuxObsNP30
     */
    double computeThValue();
      
private:
    const NPbase * myNPbase;
};




//-----  Deviations of SM inputs with respect to reference value  ----------
// (To use in future collider studies where the predictions are assumed to be SM at some
//  reference point)

/**
 * @class dalphaMzRef
 * @brief An observable class for the (relative) deviation of @f$\alpha(M_z)@f$ with respect to the SM reference value.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the (relative) deviation of @f$\alpha(M_z)@f$ with respect to the SM reference value.
 *
 */
class dalphaMzRef : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    dalphaMzRef(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the dalphaMzRef class.
     */
    virtual ~dalphaMzRef();

    /**
     * @brief The (relative) deviation of @f$\alpha(M_z)@f$ with respect to the SM reference value.
     * @return @f$\delta \alpha(M_z) / \alpha(M_z)@f$
     */
    double computeThValue();

    const NPbase * myNPbase;
    
private:


};


/**
 * @class dalphaSMzRef
 * @brief An observable class for the (relative) deviation of @f$\alpha_s(M_z)@f$ with respect to the SM reference value.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the (relative) deviation of @f$\alpha_s(M_z)@f$ with respect to the SM reference value.
 *
 */
class dalphaSMzRef : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    dalphaSMzRef(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the dalphaSMzRef class.
     */
    virtual ~dalphaSMzRef();

    /**
     * @brief The (relative) deviation of @f$\alpha_s(M_z)@f$ with respect to the SM reference value.
     * @return @f$\delta \alpha_s(M_z) / \alpha_s(M_z)@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:

    
};


/**
 * @class dMzRef
 * @brief An observable class for the (relative) deviation of @f$M_z@f$ with respect to the SM reference value.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the (relative) deviation of @f$M_z@f$ with respect to the SM reference value.
 *
 */
class dMzRef : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    dMzRef(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the dMzRef class.
     */
    virtual ~dMzRef();

    /**
     * @brief The (relative) deviation of @f$M_z@f$ with respect to the SM reference value.
     * @return @f$\delta M_z / M_z@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:

    
};


/**
 * @class dMHRef
 * @brief An observable class for the (relative) deviation of @f$M_H@f$ with respect to the SM reference value.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the (relative) deviation of @f$M_H@f$ with respect to the SM reference value.
 *
 */
class dMHRef : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    dMHRef(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the dMHRef class.
     */
    virtual ~dMHRef();

    /**
     * @brief The (relative) deviation of @f$M_H@f$ with respect to the SM reference value.
     * @return @f$\delta M_H / M_H@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:

    
};


/**
 * @class dmtRef
 * @brief An observable class for the (relative) deviation of @f$m_t@f$ with respect to the SM reference value.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the (relative) deviation of @f$m_t@f$ with respect to the SM reference value.
 *
 */
class dmtRef : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    dmtRef(const StandardModel& SM_i);
      
    /**
     * @brief Destructor of the dmtRef class.
     */
    virtual ~dmtRef();

    /**
     * @brief The (relative) deviation of @f$m_t@f$ with respect to the SM reference value.
     * @return @f$\delta m_t / m_t@f$
     */
    double computeThValue();
      
    const NPbase * myNPbase;
    
private:

    
};

// Top Wilson coefficients in the notation of LHC Top WG arXiv: 1802.07237

/**
 * @class TWGcQQ1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQQ1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQQ1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQQ1 class.
     */
    virtual ~TWGcQQ1();

    /**
     * @brief The Wilson coefficient @f$C_{QQ}^{1}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{QQ}^{1}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQQ8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQQ8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQQ8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQQ8 class.
     */
    virtual ~TWGcQQ8();

    /**
     * @brief The Wilson coefficient @f$C_{QQ}^{8}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{QQ}^{8}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctt1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctt1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctt1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctt1 class.
     */
    virtual ~TWGctt1();

    /**
     * @brief The Wilson coefficient @f$C_{tt}^{1}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{tt}^{1}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQt1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQt1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQt1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQt1 class.
     */
    virtual ~TWGcQt1();

    /**
     * @brief The Wilson coefficient @f$C_{Qt}^{1}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Qt}^{1}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQt8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQt8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQt8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQt8 class.
     */
    virtual ~TWGcQt8();

    /**
     * @brief The Wilson coefficient @f$C_{Qt}^{8}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Qt}^{8}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQq31
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQq31 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQq31(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQq31 class.
     */
    virtual ~TWGcQq31();

    /**
     * @brief The Wilson coefficient @f$C_{Q q}^{(3,1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q q}^{(3,1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQq38
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQq38 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQq38(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQq38 class.
     */
    virtual ~TWGcQq38();

    /**
     * @brief The Wilson coefficient @f$C_{Q q}^{(3,8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q q}^{(3,8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQq11
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQq11 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQq11(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQq11 class.
     */
    virtual ~TWGcQq11();

    /**
     * @brief The Wilson coefficient @f$C_{Q q}^{(1,1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q q}^{(1,1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQq18
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQq18 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQq18(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQq18 class.
     */
    virtual ~TWGcQq18();

    /**
     * @brief The Wilson coefficient @f$C_{Q q}^{(1,8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q q}^{(1,8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQu1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQu1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQu1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQu1 class.
     */
    virtual ~TWGcQu1();

    /**
     * @brief The Wilson coefficient @f$C_{Q u}^{(1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q u}^{(1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQu8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQu8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQu8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQu8 class.
     */
    virtual ~TWGcQu8();

    /**
     * @brief The Wilson coefficient @f$C_{Q u}^{(8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q u}^{(8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQd1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQd1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQd1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQd1 class.
     */
    virtual ~TWGcQd1();

    /**
     * @brief The Wilson coefficient @f$C_{Q d}^{(1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q d}^{(1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQd8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQd8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQd8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQd8 class.
     */
    virtual ~TWGcQd8();

    /**
     * @brief The Wilson coefficient @f$C_{Q d}^{(8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q d}^{(8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctq1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctq1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctq1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctq1 class.
     */
    virtual ~TWGctq1();

    /**
     * @brief The Wilson coefficient @f$C_{t q}^{(1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t q}^{(1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctq8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctq8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctq8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctq8 class.
     */
    virtual ~TWGctq8();

    /**
     * @brief The Wilson coefficient @f$C_{t q}^{(8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t q}^{(8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctu1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctu1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctu1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctu1 class.
     */
    virtual ~TWGctu1();

    /**
     * @brief The Wilson coefficient @f$C_{t u}^{(1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t u}^{(1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctu8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctu8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctu8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctu8 class.
     */
    virtual ~TWGctu8();

    /**
     * @brief The Wilson coefficient @f$C_{tu}^{(8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{tu}^{(8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctd1
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctd1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctd1(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctd1 class.
     */
    virtual ~TWGctd1();

    /**
     * @brief The Wilson coefficient @f$C_{t d}^{(1)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t d}^{(1)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctd8
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctd8 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctd8(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctd8 class.
     */
    virtual ~TWGctd8();

    /**
     * @brief The Wilson coefficient @f$C_{td}^{(8)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{td}^{(8)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctH
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctH(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctH class.
     */
    virtual ~TWGctH();

    /**
     * @brief The Wilson coefficient @f$C_{t phi}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t phi}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};

/**
 * @class TWGcHQm
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHQm : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHQm(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHQm class.
     */
    virtual ~TWGcHQm();

    /**
     * @brief The Wilson coefficient @f$C_{phi Q}^{-}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{phi Q}^{-}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcHQp
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHQp : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHQp(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHQp class.
     */
    virtual ~TWGcHQp();

    /**
     * @brief The Wilson coefficient @f$C_{phi Q}^{+}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{phi Q}^{+}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcHQ3
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHQ3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHQ3(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHQ3 class.
     */
    virtual ~TWGcHQ3();

    /**
     * @brief The Wilson coefficient @f$C_{phi Q}^{(3)}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{phi Q}^{(3)}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcHt
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHt : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHt(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHt class.
     */
    virtual ~TWGcHt();

    /**
     * @brief The Wilson coefficient @f$C_{phi t}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{phi t}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcHb
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHb(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHb class.
     */
    virtual ~TWGcHb();

    /**
     * @brief The Wilson coefficient @f$C_{phi b}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{phi b}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcHtb
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcHtb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcHtb(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcHtb class.
     */
    virtual ~TWGcHtb();

    /**
     * @brief The Wilson coefficient @f$Re[C_{phi tb}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{phi tb}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctW
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctW(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctW class.
     */
    virtual ~TWGctW();

    /**
     * @brief The Wilson coefficient @f$Re[C_{t W}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{t W}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGImctW
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGImctW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGImctW(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGImctW class.
     */
    virtual ~TWGImctW();

    /**
     * @brief The Wilson coefficient @f$Im[C_{t W}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Im[C_{t W}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctZ
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctZ(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctZ class.
     */
    virtual ~TWGctZ();

    /**
     * @brief The Wilson coefficient @f$Re[C_{t Z}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{t Z}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGImctZ
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGImctZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGImctZ(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGImctZ class.
     */
    virtual ~TWGImctZ();

    /**
     * @brief The Wilson coefficient @f$Im[C_{t Z}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Im[C_{t Z}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctG
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctG : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctG(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctG class.
     */
    virtual ~TWGctG();

    /**
     * @brief The Wilson coefficient @f$Re[C_{t G}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{t G}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcbW
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcbW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcbW(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcbW class.
     */
    virtual ~TWGcbW();

    /**
     * @brief The Wilson coefficient @f$Re[C_{b W}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{b W}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQlM
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQlM : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQlM(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQlM class.
     */
    virtual ~TWGcQlM();

    /**
     * @brief The Wilson coefficient @f$C_{Q l}^{-}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q l}^{-}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQlP
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQlP : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQlP(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQlP class.
     */
    virtual ~TWGcQlP();

    /**
     * @brief The Wilson coefficient @f$C_{Q l}^{+}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q l}^{+}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQl3
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQl3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQl3(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQl3 class.
     */
    virtual ~TWGcQl3();

    /**
     * @brief The Wilson coefficient @f$C_{Ql}^{3}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Ql}^{3}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcQe
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcQe : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcQe(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcQe class.
     */
    virtual ~TWGcQe();

    /**
     * @brief The Wilson coefficient @f$C_{Q e}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{Q e}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctl
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctl : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctl(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctl class.
     */
    virtual ~TWGctl();

    /**
     * @brief The Wilson coefficient @f$C_{t l}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t l}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGcte
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGcte : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGcte(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGcte class.
     */
    virtual ~TWGcte();

    /**
     * @brief The Wilson coefficient @f$C_{t e}/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$C_{t e}/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


/**
 * @class TWGctlS
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctlS : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctlS(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctlS class.
     */
    virtual ~TWGctlS();

    /**
     * @brief The Wilson coefficient @f$Re[C_{t}^{S(l)}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{t}^{S(l)}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};

/**
 * @class TWGctlT
 * @brief An observable class for the LHC Top WG definitions of Top operators.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the Wilson coefficient of top operators in the notation of LHC Top WG arXiv: 1802.07237.
 *
 */
class TWGctlT : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    TWGctlT(const StandardModel& SM_i, const double mu_i);
      
    /**
     * @brief Destructor of the TWGctlT class.
     */
    virtual ~TWGctlT();

    /**
     * @brief The Wilson coefficient @f$Re[C_{t}^{T(l)}]/\Lambda^2@f$ in units of @f$TeV^{-2}@f$.
     * @return @f$Re[C_{t}^{T(l)}]/\Lambda^2@f$
     */
    double computeThValue();
          
private:
    const NPbase * myNPbase;
    const double mu;

};


#endif	/* NPCOUPLINGS_H */

