#pragma once

#include <string_view>
#include <vector>

namespace optics {

struct ClusterPublisher {
    virtual ~ClusterPublisher() = 0;

    virtual void publish(std::vector<char const *> const &cluster) = 0;
};

}  // namespace optics
