/*
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSTHOBSERVABLES_H
#define	HIGGSTHOBSERVABLES_H

#include <ThObservable.h>

/**
 * @class BrWWRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f$BR(H\to WW@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$BR(H\to WW@f$
 * in the current model and in the Standard Model
 */
class BrWWRatio : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    BrWWRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrWWRatio called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the the ratio of the @f$BR(H\to WW@f$ in the current model and SM
     * @return
     */
    double computeThValue()
    {
        return myNPbase->computeKW() * myNPbase->computeKW() / myNPbase->computeGTotalRatio();
    }

private:
    const NPbase* myNPbase;
}; // either this or put back computeKW() and friends in the StandardModel class!!!

/**
 * @class BrZZRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f$BR(H\to ZZ@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$BR(H\to ZZ@f$
 * in the current model and in the Standard Model
 */
class BrZZRatio : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrZZRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrZZRatio called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the the ratio of the @f$BR(H\to ZZ@f$ in the current model and SM
     * @return
     */
    double computeThValue()
    {
        return myNPbase->computeKZ() * myNPbase->computeKZ() / myNPbase->computeGTotalRatio();
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class BrgagaRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f$BR(H\to\gamma\gamma@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$BR(H\to\gamma\gamma@f$
 * in the current model and in the Standard Model
 */
class BrgagaRatio : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    BrgagaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrgagaRatio called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the the ratio of the @f$BR(H\to\gamma\gamma@f$ in the current model and SM
     * @return
     */
    double computeThValue()
    {
        return myNPbase->computeKgaga() * myNPbase->computeKgaga() / myNPbase->computeGTotalRatio();
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class BrtautauRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the @f$BR(H\to\tau\tau@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the @f$BR(H\to\tau\tau@f$
 * in the current model and in the Standard Model
 */
class BrtautauRatio : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    BrtautauRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrtautauRatio called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the the ratio of the @f$BR(H\to\tau\tau@f$ in the current model and SM
     * @return
     */
    double computeThValue()
    {
        return myNPbase->computeKtau() * myNPbase->computeKtau() / myNPbase->computeGTotalRatio();
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muVBF
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{VBF}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */
class muVBF : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muVBF(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBF called with a class whose parent is not NPbase");

    }

    /**
     * method to compute the value of @f$\mu_{VBF}@f$ in the current model
     * @return 
     */
    double computeThValue()
    {
        return (myNPbase->computeKW() * myNPbase->computeKW() * myNPbase->computeSigmaWF() + myNPbase->computeKZ() * myNPbase->computeKZ() * myNPbase->computeSigmaZF() +
                myNPbase->computeKW() * myNPbase->computeKZ() * myNPbase->computeSigmaZWF()) /
                (myNPbase->computeSigmaWF() + myNPbase->computeSigmaZF() + myNPbase->computeSigmaZWF());
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muWH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{WH}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{WH}@f$ between the W Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muWH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muWH(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muWH called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the value of @f$\mu_{WH}@f$ in the current model
     * @return 
     */
    double computeThValue()
    {
        return (myNPbase->computeKW() * myNPbase->computeKW());
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muZH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ZH}@f$
 * @author SusyFit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ZH}@f$ between the Z Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muZH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muZH(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muZH called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the value of @f$\mu_{ZH}@f$ in the current model
     * @return 
     */
    double computeThValue()
    {
        return (myNPbase->computeKZ() * myNPbase->computeKZ());
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muVH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{VH}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VH}@f$ between the WH+ZH associated production cross-section
 * in the current model and in the Standard Model
 */
class muVH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muVH(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVH called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the value of @f$\mu_{VH}@f$ in the current model
     * @return
     */
    double computeThValue()
    {
        if (myNPbase->computeKW() == myNPbase->computeKZ())
            return (myNPbase->computeKW() * myNPbase->computeKW());
        else {
            double sigmaWH_SM = myNPbase->getTrueSM().computeSigmaWH();
            double sigmaZH_SM = myNPbase->getTrueSM().computeSigmaZH();
            return ((myNPbase->computeKW() * myNPbase->computeKW() * sigmaWH_SM
                    + myNPbase->computeKZ() * myNPbase->computeKZ() * sigmaZH_SM)
                    / (sigmaWH_SM + sigmaZH_SM));
        }
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muggH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ggH}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs production cross-section
 * in the current model and in the Standard Model
 */
class muggH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muggH(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggH called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the value of @f$\mu_{ggH}@f$ in the current model
     * @return 
     */
    double computeThValue()
    {
        return myNPbase->computeKglgl() * myNPbase->computeKglgl();
    }
private:
    const NPbase* myNPbase;
};

/**
 * @class muttH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ttH}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated production cross-section
 * in the current model and in the Standard Model
 */
class muttH : public ThObservable {
public:

    /**
     * @brief constructor
     * @param SM_i a reference to a StandardModel object or to any extension of it
     */
    muttH(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttH called with a class whose parent is not NPbase");
    }

    /**
     * method to compute the value of @f$\mu_{ttH}@f$ in the current model
     * @return 
     */
    double computeThValue()
    {
        return (myNPbase->computeKt() * myNPbase->computeKt());
    }
private:
    const NPbase* myNPbase;
};

#endif	/* HIGGSTHOBSERVABLES_H */

