/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef STD_MAKE_VECTOR_H
#define STD_MAKE_VECTOR_H

#include <vector>

template <typename T>
class make_vector {
public:
    typedef make_vector<T> my_type;

    my_type& operator<<(const T& val) {
        data_.push_back(val);
        return *this;
    }

    operator std::vector<T>() const {
        return data_;
    }
private:
    std::vector<T> data_;
};


#endif /* STD_MAKE_VECTOR_H */

