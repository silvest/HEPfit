/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSTHOBSERVABLES_H
#define	HIGGSTHOBSERVABLES_H

#include <ThObservable.h>
#include <HiggsExtensionModel.h>

/**
 * @class HiggsBaseClass
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to WW@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to WW@f]
 * in the current model and in the Standard Model
 */
class HiggsBaseClass : public ThObservable {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a const StandardModel object or to any extension of it
     */
    HiggsBaseClass(const StandardModel& HESM_i) : ThObservable(HESM_i), HESM(dynamic_cast<const HiggsExtensionModel&>(HESM_i))
    {
        if(HESM_i.ModelName().compare(0,5,"Higgs")!=0)
            throw std::runtime_error("ERROR: the HiggsBaseClass constructor can only be used with a HiggsExtensionModel reference, while I got " +
                    HESM_i.ModelName() + " as argument");
    }
        
    virtual double computeThValue()=0;
    
protected:
    const HiggsExtensionModel& HESM;
};

/**
 * @class BrWW
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to WW@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to WW@f]
 * in the current model and in the Standard Model
 */
class BrWW : public HiggsBaseClass {
public:
    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    BrWW(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}

    /**
     * method to compute the the ratio of the @f[BR(H\to WW@f] in the current model and SM
     * @return
     */    
    double computeThValue() {
        return HESM.computeKW()*HESM.computeKW()/HESM.computeGTotalRatio();
    }
    
};

/**
 * @class BrZZ
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to ZZ@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to ZZ@f]
 * in the current model and in the Standard Model
 */
class BrZZ : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrZZ(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}

    /**
     * method to compute the the ratio of the @f[BR(H\to ZZ@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKZ() * HESM.computeKZ() / HESM.computeGTotalRatio();
    }
};

/**
 * @class Brgaga
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to\gamma\gamma@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to\gamma\gamma@f]
 * in the current model and in the Standard Model
 */
class Brgaga : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    Brgaga(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}
 
    /**
     * method to compute the the ratio of the @f[BR(H\to\gamma\gamma@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKgaga()*HESM.computeKgaga()/HESM.computeGTotalRatio();
    }
};

/**
 * @class Brtautau
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f[BR(H\to\tau\tau@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f[BR(H\to\tau\tau@f]
 * in the current model and in the Standard Model
 */
class Brtautau : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    Brtautau(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}
    
    /**
     * method to compute the the ratio of the @f[BR(H\to\tau\tau@f] in the current model and SM
     * @return
     */
    double computeThValue() {
        return HESM.computeKtau()*HESM.computeKtau()/HESM.computeGTotalRatio();
    }
};

/**
 * @class muVBF
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{VBF}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{VBF}@f] between the vector-boson fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */
class muVBF : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    muVBF(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}
    
    /**
     * method to compute the value of  @f[\mu_{VBF}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKW()*HESM.computeKW()*HESM.computeSigmaWF()
                + HESM.computeKZ()*HESM.computeKZ()*HESM.computeSigmaZF()
                + HESM.computeKW()*HESM.computeKZ()*HESM.computeSigmaZWF())/
                (HESM.computeSigmaWF()+HESM.computeSigmaZF()+HESM.computeSigmaZWF());
    }
};

/**
 * @class muWH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{WH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{WH}@f] between the W Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muWH : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    muWH(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}
 
    /**
     * method to compute the value of  @f[\mu_{WH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKW()*HESM.computeKW());
    }
};

/**
 * @class muZH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ZH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ZH}@f] between the Z Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muZH : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    muZH(const HiggsExtensionModel& HESM_i) : HiggsBaseClass(HESM_i) {}
    
    /**
     * method to compute the value of  @f[\mu_{ZH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKZ()*HESM.computeKZ());
    }
};

/**
 * @class muggH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ggH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ggH}@f] between the gluon-gluon fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */
class muggH : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    muggH(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}

    /**
     * method to compute the value of  @f[\mu_{ggH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return HESM.computeKglgl()*HESM.computeKglgl();
    }
};

/**
 * @class muttH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f[\mu_{ttH}@f]
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f[\mu_{ttH}@f] between the t-tbar-Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muttH : public HiggsBaseClass {
public:

    /**
     * @brief constructor
     * @param HESM_i a reference to a StandardModel object or to any extension of it
     */
    muttH(const StandardModel& HESM_i) : HiggsBaseClass(HESM_i) {}

    /**
     * method to compute the value of  @f[\mu_{ttH}@f] in the current model
     * @return 
     */
    double computeThValue() {
        return (HESM.computeKt()*HESM.computeKt());
    }
};

#endif	/* HIGGSTHOBSERVABLES_H */

