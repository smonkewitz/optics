#include "Tree.h"

#include <absl/log/log.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>

namespace optics {

namespace {

// Finds the dimension in which the given points have maximum extent. Used to pick a
// splitting dimension during tree construction. Assumes that points != nullptr and
// numPoints > 0.
//
// Returns a pair containing the maximum extent of the input points and the dimension of
// maximum extent.
std::pair<double, size_t> MaxExtentAndDim(Point const* points, size_t numPoints) {
    Vec3 min = {std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity()};
    Vec3 max = -min;

    for (size_t i = 0; i < numPoints; ++i) {
        min = Min(min, points[i].v);
        max = Max(max, points[i].v);
    }
    Vec3 extents = max - min;
    double maxExtent = extents.x();
    size_t maxDim = 0;
    if (extents.y() > maxExtent) {
        maxExtent = extents.y();
        maxDim = 1;
    }
    if (extents.z() > maxExtent) {
        maxExtent = extents.z();
        maxDim = 2;
    }
    return std::make_pair(maxExtent, maxDim);
}

// Orders N-dimensional points along a single dimension.
struct PointCmp {
    size_t dim;

    bool operator()(Point const& a, Point const& b) const {
        return a.v.coords[dim] < b.v.coords[dim];
    }
};

}  // namespace

Tree::Tree(Point* points, size_t numPoints, size_t pointsPerLeaf,
           double leafExtentThreshold)
    : points_(points), numPoints_(numPoints) {
    if (points == nullptr || numPoints == 0) {
        throw std::invalid_argument("no input points provided");
    }
    if (pointsPerLeaf == 0) {
        throw std::invalid_argument("target number of points per leaf must be > 0");
    }
    // compute tree height
    size_t h = 0;
    while (h < MAX_HEIGHT &&
           numPoints / (static_cast<size_t>(1) << h) > pointsPerLeaf) {
        ++h;
    }
    height_ = h;
    size_t numNodes = (static_cast<size_t>(1) << (h + 1)) - 1;
    nodes_ = std::make_unique<Node[]>(numNodes);
    build(leafExtentThreshold);
}

size_t Tree::inRange(Vec3 const& v, double const dist) {
    std::array<bool, MAX_HEIGHT> descend;
    size_t node = 0;
    size_t h = 0;
    size_t head = NOT_FOUND;
    size_t tail = NOT_FOUND;
    while (true) {
        if (nodes_[node].isLeaf()) {
            // reached a leaf
            size_t left = 0;
            size_t right = nodes_[node].right();
            if ((node & (node + 1)) != 0) {
                // node has a left sibling - use it to obtain index of first point in
                // leaf.
                left = nodes_[node - 1].right();
            }
            // Scan leaf for results, and append them to embedded linked list
            for (size_t i = left; i < right; ++i) {
                double d = SquaredEuclidianDistance(v, points_[i].v);
                if (d <= dist) {
                    points_[i].dist = d;
                    if (tail == NOT_FOUND) {
                        head = i;
                    } else {
                        points_[tail].next = i;
                    }
                    tail = i;
                }
            }
            // move back up the tree
            node = (node - 1) >> 1;
            --h;
            for (; h != static_cast<size_t>(-1) && !descend[h]; --h) {
                node = (node - 1) >> 1;
            }
            if (h == static_cast<size_t>(-1)) {
                // finished tree traversal
                break;
            }
            descend[h] = false;
            node = (node << 1) + 2;
            ++h;
        } else {
            // determine which children must be visited
            double split = nodes_[node].split;
            size_t dim = nodes_[node].splitDim();
            double vd = v.coords[dim];
            if (MinSquaredEuclidianDistance(vd, split) <= dist) {
                // both children must be visited
                descend[h] = true;
                node = (node << 1) + 1;
                ++h;
            } else if (vd < split) {
                // visit left child
                descend[h] = false;
                node = (node << 1) + 1;
                ++h;
            } else {
                // visit right child
                descend[h] = false;
                node = (node << 1) + 2;
                ++h;
            }
        }
    }
    return head;
}

void Tree::build(double leafExtentThreshold) {
    LOG(INFO) << "building 3d tree of height " << height_ << " for " << numPoints_
              << " points";
    size_t node = 0;
    size_t left = 0;
    size_t right = numPoints_;
    size_t h = 0;
    while (true) {
        nodes_[node].setRight(right);
        if (h < height_) {
            // find splitting dimension
            auto [extent, dim] = MaxExtentAndDim(points_ + left, right - left);
            if (extent > leafExtentThreshold) {
                nodes_[node].setSplitDim(dim);
                // find median of array
                size_t median = left + ((right - left) >> 1);
                std::nth_element(points_ + left, points_ + median, points_ + right,
                                 PointCmp{dim});
                right = median;
                nodes_[node].split = points_[right].v.coords[dim];
                // process left child
                node = (node << 1) + 1;
                ++h;
                continue;
            }
            // node extent is below the subdivision limit: set right index for all right
            // children of node as their left siblings may be valid
            size_t h2 = h;
            size_t c = node;
            do {
                c = (c << 1) + 2;
                ++h2;
                nodes_[c].setRight(right);
            } while (h2 < height_);
        }
        // move up the tree until a left child is found
        left = right;
        for (; h > 0 && (node & 1) == 0; --h) {
            node = (node - 1) >> 1;
        }
        if (h == 0) {
            // tree construction complete!
            break;
        }
        // node is now the index of a left child - process its right sibling
        right = nodes_[(node - 1) >> 1].right();
        node += 1;
    }
    LOG(INFO) << "built 3d tree";
}

}  // namespace optics
