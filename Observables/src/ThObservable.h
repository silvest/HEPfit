/* 
 * File:   ThObservable.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:15 PM
 */

#ifndef THOBSERVABLE_H
#define	THOBSERVABLE_H

class ThObservable {
public:
    ThObservable();
    ThObservable(const ThObservable& orig);
    virtual ~ThObservable();
    virtual double getThValue() = 0;
private:
};

#endif	/* THOBSERVABLE_H */

