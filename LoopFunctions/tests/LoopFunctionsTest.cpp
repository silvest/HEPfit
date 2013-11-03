/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include "PVfunctions.h"

using namespace std;

int main(int argc, char** argv) {
    
    PVfunctions PVtest;    
    LoopToolsWrapper LT;
    
    try {
    
        double Mw2 = 80.0*80.0;
        cout << "B0p(Mw2; Mw2; Mw2,Mw2)= " << PVtest.B0p(Mw2,Mw2,Mw2,Mw2,false) << " "
             << LT.PV_B0p(Mw2,Mw2,Mw2,Mw2) << endl;
        cout << "B0p(Mw2; Mw2; 0,Mw2)= " << PVtest.B0p(Mw2,Mw2,0.0,Mw2,false) << " "
             << LT.PV_B0p(Mw2,Mw2,0.0,Mw2) << endl;
        cout << "B0p(Mw2; Mw2; 0,0)= " << PVtest.B0p(Mw2,Mw2,0.0,0.0,false) << " "
             << LT.PV_B0p(Mw2,Mw2,0.0,0.0) << endl;
        cout << "B0p(Mw2; 0; Mw2,Mw2)= " << PVtest.B0p(Mw2,0.0,Mw2,Mw2,false) << " "
             << LT.PV_B0p(Mw2,0.0,Mw2,Mw2) << endl;

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

