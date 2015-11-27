/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSTHOBSERVABLES_H
#define	HIGGSTHOBSERVABLES_H

#include "ThObservable.h"
#include "NPbase.h"

/**
 * @class muggH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ggH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ggH}@f$ between the gluon-gluon
 * fusion Higgs production cross-section in the current model and in the Standard Model.
 */
class muggH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ggH}@f$ in the current model.
     * @return @f$\mu_{ggH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muggH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBF
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{VBF}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF}@f$ between the vector-boson
 * fusion Higgs production cross-section in the current model and in the Standard Model.
 */
class muVBF : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBF(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBF called with a class whose parent is not NPbase");

    }

    /**
     * @brief A method to compute the value of @f$\mu_{VBF}@f$ in the current model.
     * @return @f$\mu_{VBF}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVBF(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muWH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{WH}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{WH}@f$ between the W Higgs 
 * associated production cross-section in the current model and in the Standard Model. 
 */
class muWH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muWH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muWH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{WH}@f$ in the current model.
     * @return @f$\mu_{WH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muWH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muZH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ZH}@f$ between the Z Higgs 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muZH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muZH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muZH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ZH}@f$ in the current model.
     * @return @f$\mu_{ZH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muZH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeZH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to ZH}@f$ between the 
 * @f$e^+e^- \to ZH@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeZH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeZH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeZH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to ZH}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to ZH}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeZH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{VH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VH}@f$ between the WH+ZH 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muVH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{VH}@f$ in the current model.
     * @return @f$\mu_{VH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muVBFpVH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{VBF+VH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{VBF+VH}@f$ between the sum of
 * VBF and WH+ZH 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muVBFpVH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muVBFpVH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muVBFpVH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{VBF+VH}@f$ in the current model.
     * @return @f$\mu_{VBF+VH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muVBFpVH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muttH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ttH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs 
 * associated production cross-section in the current model and in the Standard Model.
 */
class muttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muttH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muttH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ttH}@f$ in the current model. 
     * @return @f$\mu_{ttH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muttH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class muggHpttH
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio @f$\mu_{ggH+ttH}@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{ggH+ttH}@f$ between the sum
 * of gluon-gluon fusion and t-tbar-Higgs associated production cross-section in 
 * the current model and in the Standard Model.
 */
class muggHpttH : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    muggHpttH(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("muggHpttH called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{ggH+ttH}@f$ in the current model. 
     * @return @f$\mu_{ggH+ttH}@f$
     */
    double computeThValue()
    {
        return myNPbase->muggHpttH(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class BrHtoggRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to gg)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to gg)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoggRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoggRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoggRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to gg)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to gg)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHggRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoWWRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to WW)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to WW)@f$.
 * in the current model and in the Standard Model
 */
class BrHtoWWRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoWWRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoWWRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to WW)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to WW)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHWWRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZZRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to ZZ)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to ZZ)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZZRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a HiggsExtensionModel object or to any extension of it
     */
    BrHtoZZRatio(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZZRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to ZZ)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZZRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoZgaRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to Z\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtoZgaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoZgaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoZgaRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to Z\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHZgaRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtogagaRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to\gamma\gamma)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to\gamma\gamma)@f$
 * in the current model and in the Standard Model.
 */
class BrHtogagaRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtogagaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtogagaRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHgagaRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtotautauRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to \tau^+\tau^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to \tau^+\tau^-)@f$
 * in the current model and in the Standard Model.
 */
class BrHtotautauRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtotautauRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtotautauRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to \tau^+\tau^-)@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHtautauRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtoccRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to c\bar{c})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to c\bar{c})@f$
 * in the current model and in the Standard Model.
 */
class BrHtoccRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtoccRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtoccRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to c\bar{c})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHccRatio();
    }

private:
    const NPbase* myNPbase;
};

/**
 * @class BrHtobbRatio
 * @ingroup HiggsExtensions
 * @brief A class for computing the ratio of the Br@f$(H\to b\bar{b})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the ratio of the Br@f$(H\to b\bar{b})@f$
 * in the current model and in the Standard Model.
 */
class BrHtobbRatio : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     */
    BrHtobbRatio(const StandardModel& SM_i) : ThObservable(SM_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("BrHtobbRatio called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the the ratio of the Br@f$(H\to b\bar{b})@f$
     * in the current model and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$
     */
    double computeThValue()
    {
        return myNPbase->BrHbbRatio();
    }

private:
    const NPbase* myNPbase;
};

#endif	/* HIGGSTHOBSERVABLES_H */

