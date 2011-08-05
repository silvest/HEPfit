/* 
 * File:   LoopFunctionsTest.cpp
 * Author: mishima
 *
 * Created on Aug 5, 2011, 2:00:06 AM
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "PVfunctions.h"

using namespace std;

int main(int argc, char** argv) {
    
    PVfunctions PVtest;    
    
    try {
    
        // Invalid
        //cout << PVtest.B0(-1.0,1.0,1.0,1.0) << endl;
        //cout << PVtest.B0(1.0,0.0,0.0,0.0) << endl;

        cout << "B0(5,1,2,3)= " << PVtest.B0(5.0,1.0,2.0,3.0) << endl;
        cout << "B0(1,1,1,1)= " << PVtest.B0(1.0,1.0,1.0,1.0) << endl;
        cout << "B0(1,0,1,1)= " << PVtest.B0(1.0,0.0,1.0,1.0) << endl;        
        cout << "B0(1,1,0,1)= " << PVtest.B0(1.0,1.0,0.0,1.0) << endl;
        cout << "B0(1,2,0,1)= " << PVtest.B0(1.0,2.0,0.0,1.0) << endl;
        cout << "B0(1,1,1,0)= " << PVtest.B0(1.0,1.0,1.0,0.0) << endl;
        cout << "B0(1,0,0,1)= " << PVtest.B0(1.0,0.0,0.0,1.0) << endl;
        cout << "B0(1,1,0,0)= " << PVtest.B0(1.0,1.0,0.0,0.0) << endl;
    
        
        return EXIT_SUCCESS;
    } catch (const char* c) {
        cerr << c << endl;
        return EXIT_FAILURE;
    }    
}

