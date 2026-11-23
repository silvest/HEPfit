/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MTMSBAR_H
#define MTMSBAR_H

#include "ThObservable.h"

/**
 * @class MtMSbar
 * @ingroup StandardModel
 * @brief A class for @f$m_t(m_t)@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the value of
 * @f$m_t(m_t)@f$ in the @f$\bar{\mathrm{MS}}@f$ scheme from the top pole mass, using the three-loop
 * (FULLNNNLO) pole @f$\leftrightarrow\bar{\mathrm{MS}}@f$ relation, the same order used by QCD::PostUpdate()
 * to convert between the top input mass and the pole mass, so the value coincides with the internally used
 * @f$m_t(m_t)@f$.
 */
class MtMSbar : public ThObservable {
public:
    /**
     * @brief Constructor. 
     * @param[in] SM_i
     */
    MtMSbar(const StandardModel& SM_i) : ThObservable(SM_i) {};

    /**
     *
     * @return theoretical value of @f$\bar{\mathrm{MS}}@f$
     */
    virtual double computeThValue();

    private:

};

#endif /* MTMSBAR_H */

