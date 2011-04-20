/* 
 * File:   ThFactory.cpp
 * Author: silvest
 * 
 * Created on April 19, 2011, 12:23 PM
 */

#include "ThFactory.h"
#include <Flavour.h>

ThFactory::ThFactory(StandardModel * myModel) {
    thobs["Dmd0"] = new Dmb(myModel,0);
    thobs["Dmd1"] = new Dmb(myModel,1);
    thobs["Vud"] = new Vud(myModel);
    thobs["Vus"] = new Vus(myModel);
    thobs["Vub"] = new Vub(myModel);
    thobs["Vcb"] = new Vcb(myModel);
    thobs["alpha"] = new alpha(myModel);
    thobs["alpha_2a"] = new alpha_2a(myModel);
    thobs["gamma"] = new gammac(myModel);
    }

ThFactory::ThFactory(const ThFactory& orig) {
    thobs = orig.thobs;
}

ThFactory::~ThFactory() {
    for(std::map<std::string,ThObservable *>::iterator it = thobs.begin(); 
            it != thobs.end(); it++)
        if(it->second != NULL) 
            delete it->second;
}

ThObservable * ThFactory::getThMethod(const std::string& name){
    if(thobs.find(name) == thobs.end()){
        std::cout << "wrong observable " << name <<" in ThFactory" << std::endl;
        exit(EXIT_FAILURE);        
    }
    return(thobs[name]);
}
