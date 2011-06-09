/* 
 * File:   Abottom.cpp
 * Author: mishima
 * 
 * Created on June 9, 2011, 3:43 PM
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Abottom.h"

Abottom::Abottom(const EW& myEW) : ThObservable(myEW) {
}

double Abottom::getThValue() {   
    
    std::cout << "Write codes!" << std::endl;
    exit(EXIT_FAILURE); 
    
    return (0.0);
}
        

