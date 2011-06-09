/* 
 * File:   sigmaHadron.cpp
 * Author: mishima
 * 
 * Created on June 9, 2011, 3:44 PM
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "sigmaHadron.h"

sigmaHadron::sigmaHadron(const EW& myEW) : ThObservable(myEW) {
}

double sigmaHadron::getThValue() {   
    
    std::cout << "Write codes!" << std::endl;
    exit(EXIT_FAILURE); 
    
    return (0.0);
}
        


