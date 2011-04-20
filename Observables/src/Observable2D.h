/* 
 * File:   Observable2D.h
 * Author: silvest
 *
 * Created on April 19, 2011, 3:36 PM
 */

#ifndef OBSERVABLE2D_H
#define	OBSERVABLE2D_H

#include "Observable.h"

class Observable2D : public Observable {
public:
    Observable2D(const std::string name_i, const std::string thname_i,
        const std::string thname2_i, const std::string label_i, 
        const std::string label2_i, const bool tMCMC_i, const double min_i,
        const double max_i, const double min2_i, const double max2_i, 
        ThObservable * tho_i, ThObservable * tho2_i);
    Observable2D(const Observable& o1d);
    Observable2D(const Observable2D& orig);
    virtual ~Observable2D();

    double getTheoryValue2();

    std::string getLabel2() const {
        return label2;
    }

    void setLabel2(std::string label2) {
        this->label2 = label2;
    }

    double getMax2() const {
        return max2;
    }

    void setMax2(double max2) {
        this->max2 = max2;
    }

    double getMin2() const {
        return min2;
    }

    void setMin2(double min2) {
        this->min2 = min2;
    }

    std::string getThname2() const {
        return thname2;
    }

    void setThname2(std::string thname2) {
        this->thname2 = thname2;
    }

    ThObservable* getTho2() const {
        return tho2;
    }

    void setTho2(ThObservable* tho2) {
        this->tho2 = tho2;
    }

private:
    std::string thname2, label2;
    double min2, max2;
    ThObservable * tho2;
};

#endif	/* OBSERVABLE2D_H */

