/* 
 * File:   AFBlepton.cpp
 * Author: mishima
 * 
 * Created on June 9, 2011, 3:41 PM
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "AFBlepton.h"

AFBlepton::AFBlepton(const EW& myEW) : ThObservable(myEW) {
}

double AFBlepton::getThValue() {   
    
    std::cout << "Write codes!" << std::endl;
    exit(EXIT_FAILURE); 
    
    return (0.0);
}
        

