/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMMATCHING_H
#define	GMMATCHING_H

#include "StandardModelMatching.h"

class GeorgiMachacek;

/**
 * @class GMMatching
 * @ingroup GeorgiMachacek
 * @brief A class for the Wilson coefficients in the Georgi-Machacek model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is empty.
 */
class GMMatching : public StandardModelMatching {
public:
    GMMatching(const GeorgiMachacek & GeorgiMachacek_i);

private:
    const GeorgiMachacek & myGM;
};

#endif	/* GMMATCHING_H */
