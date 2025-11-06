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

class Vcb : public ThObservable
{
    public:

        /**
         * @brief Constructor. 
         * @details This isn't needed in this detail, but will be useful later on with the generalization
         * @param[in] SM_i a reference to an object of type StandardModel
         * @param[in] meson_i initial meson of the decay
         */

        Vcb(const StandardModel& SM_i, QCD::meson meson_i);

        /**
         * @brief The computation of the moments, varying the number and elcut
         * @note All the other parameters are initialized in the constructor
         */

        double computeThValue();

    private:
        QCD::meson meson1; /**< Initial meson type */

};


#endif /* BCLNUOBSERVABLES_H */