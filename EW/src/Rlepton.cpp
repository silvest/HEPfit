/* 
 * File:   Rlepton.cpp
 * Author: mishima
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Rlepton.h"

Rlepton::Rlepton(const EW& EW_i) : ThObservable(EW_i) {
}

double Rlepton::getThValue() {   
    
    std::cout << "Write codes!" << std::endl;
    exit(EXIT_FAILURE); 
    
    return (0.0);
}
        

