/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelParameter.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

ModelParameter::ModelParameter(std::string name_in, double ave_in, double errg_in, double errf_in) 
{
    name = name_in;
    ave = ave_in;
    errg = errg_in;
    errf = errf_in;
    min = ave - errf - 5.*errg;
    max = ave + errf + 5.*errg;
    cgp_name.clear();
    if(errg == 0. && errf == 0.)
        isFixed = true;
    else
        isFixed = false;
}

ModelParameter::ModelParameter(std::string name_in, double central, Intervals intervals)
    : name(std::move(name_in)), ave(central), errg(0.), errf(0.), isFixed(false),
      flatIntervals(std::move(intervals))
{
    const std::string context = "ModelParameter " + name + " flat_intervals: ";
    if (flatIntervals.size() < 2)
        throw std::runtime_error(context + "expected at least two intervals");
    bool containsCentral = false;
    for (const auto& interval : flatIntervals) {
        if (!std::isfinite(interval.first) || !std::isfinite(interval.second) || interval.first >= interval.second)
            throw std::runtime_error(context + "each interval needs finite endpoints with lower < upper");
        if (interval.first <= central && central <= interval.second)
            containsCentral = true;
    }
    if (!containsCentral)
        throw std::runtime_error(context + "central value must belong to the support");
    std::sort(flatIntervals.begin(), flatIntervals.end());
    for (size_t i = 1; i < flatIntervals.size(); ++i)
        if (flatIntervals[i].first <= flatIntervals[i - 1].second)
            throw std::runtime_error(context + "intervals must be separated by a gap");
    min = flatIntervals.front().first;
    max = flatIntervals.back().second;
    if (!std::isfinite(max - min))
        throw std::runtime_error(context + "parameter range is too wide");
}

ModelParameter::ModelParameter() 
{}

ModelParameter::~ModelParameter() 
{
}

std::ostream& operator<<(std::ostream& output, const ModelParameter& m)
{
    if (!m.flatIntervals.empty()) {
        output << "ModelParameter " << m.name << " " << m.ave << " flat_intervals";
        for (const auto& interval : m.flatIntervals)
            output << " " << interval.first << " " << interval.second;
        return output << std::endl;
    }
    output << "ModelParameter name, average, gaussian error, flat error" << std::endl;
    output << m.name << " " << m.ave << " " << m.errg << " " << m.errf << std::endl;
    return output;
}

boost::tokenizer<boost::char_separator<char> >::iterator & ModelParameter::ParseModelParameter(boost::tokenizer<boost::char_separator<char> >::iterator & beg)
{
    auto next = [&beg]() {
        if (beg.at_end() || beg->front() == '#')
            throw std::runtime_error("Missing field in ModelParameter");
        return *beg++;
    };
    const std::string parameterName = next();
    const std::string central = next();
    const std::string prior = next();
    if (prior == "flat_intervals") {
        auto number = [&parameterName](const std::string& field) {
            try {
                return boost::lexical_cast<double>(field);
            } catch (const boost::bad_lexical_cast&) {
                throw std::runtime_error("ModelParameter " + parameterName + ": invalid number '" + field + "' in flat_intervals");
            }
        };
        Intervals intervals;
        while (!beg.at_end() && beg->front() != '#') {
            const double low = number(next());
            const double high = number(next());
            intervals.emplace_back(low, high);
        }
        *this = ModelParameter(parameterName, number(central), std::move(intervals));
    } else {
        const std::string flatError = next();
        *this = ModelParameter(parameterName, std::atof(central.c_str()),
                               std::atof(prior.c_str()), std::atof(flatError.c_str()));
    }
    return beg;
}
