#include "ModelParameter.h"

#include <iostream>
#include <stdexcept>
#include <string>

static ModelParameter parse(const std::string& fields)
{
    boost::char_separator<char> separator(" \t");
    boost::tokenizer<boost::char_separator<char> > tokens(fields, separator);
    auto it = tokens.begin();
    ModelParameter result;
    result.ParseModelParameter(it);
    return result;
}

static void check(bool condition)
{
    if (!condition) throw std::runtime_error("ModelParameter config assertion failed");
}

static void reject(const std::string& fields)
{
    try {
        parse(fields);
    } catch (const std::runtime_error&) {
        return;
    }
    throw std::runtime_error("Accepted invalid ModelParameter: " + fields);
}

int main()
{
    const ModelParameter parameter = parse("x 1 flat_intervals 10 13 -2 -1 0 1 # three components");
    const ModelParameter::Intervals expected{{-2, -1}, {0, 1}, {10, 13}};
    check(parameter.getFlatIntervals() == expected);
    check(parameter.getave() == 1 && parameter.getmin() == -2 && parameter.getmax() == 13);
    check(!parameter.IsFixed() && !parameter.IsCorrelated());
    check(!parse("x -2 flat_intervals -2 -1 3 9").IsFixed());

    reject("x 0 flat_intervals");
    reject("x 0 flat_intervals -1 1");
    reject("x 0 flat_intervals -1 1 2");
    reject("x 0 flat_intervals -1 1 2 # missing endpoint");
    reject("x 0 flat_intervals -1 1 0 2");
    reject("x 0 flat_intervals -1 1 1 2");
    reject("x 0 flat_intervals -1 1 3 2");
    reject("x 0 flat_intervals -1 1 2 2");
    reject("x 2 flat_intervals -1 1 3 4");
    reject("x nan flat_intervals -1 1 2 3");
    reject("x 0 flat_intervals -1 1 nan 3");
    reject("x 0 flat_intervals -1 1 2 inf");
    reject("x 0 flat_intervals -1 1 2 3x");
    reject("x 0 flat_intervals -1e308 0 1 1e308");

    const ModelParameter legacy = parse("x 3 1 2");
    check(legacy.getave() == 3 && legacy.geterrg() == 1 && legacy.geterrf() == 2);
    check(legacy.getmin() == -4 && legacy.getmax() == 10 && legacy.getFlatIntervals().empty());
    check(parse("x 1 0 0").IsFixed());
    std::cout << "ModelParameter flat interval syntax, support and legacy priors validated\n";
}
