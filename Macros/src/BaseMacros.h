/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <TString.h>
#include <TText.h>
#include <TLatex.h>
#include <TColor.h>

#ifndef BASEMACROS_H
#define	BASEMACROS_H

/**
 * @class BaseMacros
 * @ingroup Macros
 * @brief A base class for drawing histograms. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class BaseMacros {
public:
    
    BaseMacros()
    {
    }

    /*
     * Define new colors. 
     */
    static void DefineNewColours();

    /**
     * Convert the names of parameters and observables into the ROOT form. 
     */
    TString ConvertTitle(TString orig) const;


private:

};

#endif	/* BASEMACROS_H */

