#pragma once

#include <cstddef>
#include <memory>

#include "ClusterPublisher.h"
#include "SeedList.h"
#include "Tree.h"

namespace optics {

// Implements the OPTICS algorithm. For details, see the following paper:
//
// "OPTICS: Ordering Points To Identify the Clustering Structure".
// Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jorg Sander (1999).
// ACM SIGMOD international conference on Management of data.
// ACM Press. pp. 49-60.
class Optics {
   public:
    Optics(Point* points, size_t numPoints, size_t minNeighbors, double epsilon,
           double leafExtentThreshold, size_t pointsPerLeaf);

    void run(ClusterPublisher& publisher);

   private:
    Point* points_;  // unowned
    size_t numPoints_;
    Tree tree_;
    SeedList seeds_;
    std::unique_ptr<double[]> distances_;
    double epsilon_;
    size_t minNeighbors_;

    void expandClusterOrder(size_t i);
};

}  // namespace optics
