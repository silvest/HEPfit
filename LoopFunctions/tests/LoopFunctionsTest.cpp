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
    
    try {
    
        // Invalid
        //cout << PVtest.B0(-1.0,1.0,1.0,1.0) << endl;
        //cout << PVtest.B0(1.0,0.0,0.0,0.0) << endl;

        //cout << "B0(5,1,2,3)= " << PVtest.B0(5.0,1.0,2.0,3.0) << endl;
        //cout << "B0(1,1,1,1)= " << PVtest.B0(1.0,1.0,1.0,1.0) << endl;
        //cout << "B0(1,0,1,1)= " << PVtest.B0(1.0,0.0,1.0,1.0) << endl;        
        //cout << "B0(1,1,0,1)= " << PVtest.B0(1.0,1.0,0.0,1.0) << endl;
        //cout << "B0(1,2,0,1)= " << PVtest.B0(1.0,2.0,0.0,1.0) << endl;
        //cout << "B0(1,1,1,0)= " << PVtest.B0(1.0,1.0,1.0,0.0) << endl;
        //cout << "B0(1,0,0,1)= " << PVtest.B0(1.0,0.0,0.0,1.0) << endl;
        //cout << "B0(1,1,0,0)= " << PVtest.B0(1.0,1.0,0.0,0.0) << endl;

        double Mw = 80.36249466866697, s, sqrts;
        for (int i=0; i<100; i++) {
            sqrts = 158.0837322478437 - 0.001 + 0.002*(double)i/100.0;
            s = sqrts*sqrts;
            cout << setprecision(16) << sqrts << " " << PVtest.B0(Mw,s,Mw,Mw) << endl;

        }

        cout << PVtest.B0(Mw,158.0837322478437*158.0837322478437,Mw,Mw) << endl;
        

        //double Mz = 90.;
        //double MHp = 600., MA=900.,MH0=580.;
        //cout << PVtest.B22(Mz,Mz*Mz,MH0,MA)/Mz/Mz/M_PI << endl;
        //cout << PVtest.B22(Mz,0.0,MH0,MA)/Mz/Mz/M_PI << endl;        
        //cout << PVtest.B22(Mz,Mz*Mz,MHp,MHp)/Mz/Mz/M_PI << endl;
        //cout << PVtest.B22(Mz,0.0,MHp,MHp)/Mz/Mz/M_PI << endl;        
        
        
        //PVfunctions PVtest2;
        
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

