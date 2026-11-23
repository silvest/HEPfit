/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * (Made by Danilo Ricciardella)
 */

#ifndef BCLNUOBSERVABLES_H
#define BCLNUOBSERVABLES_H

#include "QCD.h"
#include "ThObservable.h"

class Q2moments_BClnu : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         * @param[in] typ_i the order of the moment
         */

        Q2moments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i);

        /**
         * @brief The computation of the moments, varying the number and q2cut
         * @return Q2moments(i,q2cut)
         * @note All the other parameters are initialized in the constructor
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */
        unsigned int typ; /**< Type of the moment */
};

class Elmoments_BClnu : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         * @param[in] typ_i the order of the moment
         */

        Elmoments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i);

        /**
         * @brief The computation of the moments, varying the number and elcut
         * @return Elmoments(i,elcut)
         * @note All the other parameters are initialized in the constructor
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */
        unsigned int typ; /**< Type of the moment */
};

class MXmoments_BClnu : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         * @param[in] typ_i the order of the moment
         */

        MXmoments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i);

        /**
         * @brief The computation of the moments, varying the number and elcut
         * @return Elmoments(i,elcut)
         * @note All the other parameters are initialized in the constructor
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */
        unsigned int typ; /**< Type of the moment */
};

class PartialAverageBR_BClnu : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         */

        PartialAverageBR_BClnu(const StandardModel& SM_i, QCD::meson meson_i);

        /**
         * @brief The computation of the moments, varying the number and elcut
         * @note All the other parameters are initialized in the constructor
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */
};

class Vcb_BClnu : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         */

        Vcb_BClnu(const StandardModel& SM_i, QCD::meson meson_i);

        /**
         * @brief \f$ |V_{cb}| = \sqrt{\mathcal{B}\,\hbar/(\Gamma_0 X \tau_B)} \f$ from the inclusive
         * \f$ B\to X_c\ell\nu \f$ branching ratio.
         * @details \f$ \tau_B \f$ is the average of the \f$ B^0 \f$ and \f$ B^+ \f$ lifetimes, taken from the
         * model parameters tBd and tBp (which the configuration must therefore provide, together with MBp).
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */

};


#endif /* BCLNUOBSERVABLES_H */