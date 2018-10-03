#ifndef OPT_H
#define OPT_H

#include <vector>

namespace Opt {
    enum class SearchType {SearchMinimum, SearchMaximum};

    using Point = std::vector<double>;
    using TargetFunction = double(*)(const std::vector<double>&);
}

#endif // OPT_H
