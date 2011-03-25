/* 
 * File:   Observable.h
 * Author: silvest
 *
 * Created on February 22, 2011, 11:45 AM
 */

#ifndef OBSERVABLE_H
#define	OBSERVABLE_H

#include <string>
#include <iostream>

class Observable {
public:
    Observable(const std::string name, const bool tMCMC, const double min,
    const double max, const std::string distr);
    Observable(const std::string name, const bool tMCMC, const double min,
    const double max, const std::string distr, const std::string filename);
    Observable(const std::string name, const bool tMCMC, const double min,
    const double max, const std::string distr, const double ave,
    const double errg, const double errf);
    Observable(const Observable& orig);
    virtual double getTheoryValue() = 0;
    virtual ~Observable();
//    void computeEvent(const Parameters&) const;
//    void SetHistogram(const Histopar&) const;
    std::string name, distr, filename;
    double ave,errg,errf,min,max;
    bool tMCMC;
    friend std::ostream& operator<<(std::ostream& output, const Observable& o);
private:
    void Init(const std::string name, const bool tMCMC, const double min,
    const double max, const std::string distr, const std::string filename, const
    double ave, const double errg, const double errf);
};

#endif	/* OBSERVABLE_H */

