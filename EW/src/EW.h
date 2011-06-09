/* 
 * File:   EW.h
 * Author: mishima
 *
 * Created on June 8, 2011, 4:21 PM
 */

#ifndef EW_H
#define	EW_H

#include <ThObsType.h>
#include <StandardModel.h>

using namespace gslpp;

class EW : public ThObsType {
public:

    /**
     * @brief EW constructor
     * @param[in] mySM an object of StandardModel class
     */
    EW(const StandardModel& mySM);
    
private:
    
};

#endif	/* EW_H */

