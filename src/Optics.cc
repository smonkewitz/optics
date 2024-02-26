#include "Optics.h"

#include <absl/log/check.h>
#include <absl/log/log.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace optics {

Optics::Optics(Point *points, size_t numPoints, size_t minNeighbors, double epsilon,
               double leafExtentThreshold, size_t pointsPerLeaf)
    : points_{points},
      numPoints_{numPoints},
      tree_{points, numPoints, pointsPerLeaf, leafExtentThreshold},
      seeds_{points, numPoints},
      epsilon_{std::abs(epsilon)},
      minNeighbors_{minNeighbors} {}

void Optics::run(ClusterPublisher &publisher) {
    if (points_ == nullptr) {
        throw std::runtime_error("OPTICS has already been run");
    }

    LOG(INFO) << "clustering " << numPoints_ << " points using OPTICS";
    std::vector<char const *> cluster;
    size_t scanFrom = 0;

    while (true) {
        size_t i;
        if (seeds_.empty()) {
            // find next unprocessed point
            for (i = scanFrom; i < numPoints_; ++i) {
                if (points_[i].state == UNPROCESSED) {
                    scanFrom = i + 1;
                    break;
                }
            }
            if (i == numPoints_) {
                break;
            }
            points_[i].state = PROCESSED;
            expandClusterOrder(i);
            if (cluster.size() > 0) {
                // clusters of size 1 are generated for noise sources
                publisher.publish(cluster);
                cluster.clear();
            }
            cluster.push_back(points_[i].record);
        } else {
            // expand cluster around seed with smallest reachability-distance
            i = seeds_.pop();
            expandClusterOrder(i);
            DCHECK(points_[i].reach != std::numeric_limits<double>::infinity());
            cluster.push_back(points_[i].record);
        }
    }

    publisher.publish(cluster);
    points_ = nullptr;
    LOG(INFO) << "finished clustering";
}

void Optics::expandClusterOrder(size_t i) {
    // find epsilon neighborhood of point i
    size_t const range = tree_.inRange(points_[i].v, epsilon_);
    // compute core-distance
    size_t n = 0;
    size_t j = range;
    while (j != NOT_FOUND) {
        Point *p = points_ + j;
        if (j != i) {
            double d = p->dist;
            if (n < minNeighbors_) {
                distances_[n++] = d;
                std::push_heap(distances_.get(), distances_.get() + n);
            } else if (distances_[0] > d) {
                std::pop_heap(distances_.get(), distances_.get() + n);
                distances_[n - 1] = d;
                std::push_heap(distances_.get(), distances_.get() + n);
            }
        }
        j = p->next;
    }
    if (n == minNeighbors_) {
        // point i is a core-object. Update reachability-distance of all points in the
        // epsilon-neighborhood of point i.
        double const coreDist = distances_[0];
        j = range;
        while (j != NOT_FOUND) {
            Point *p = points_ + j;
            if (p->state != PROCESSED) {
                seeds_.update(j, std::max(coreDist, p->dist));
            }
            j = p->next;
        }
    }
}

}  // namespace optics
