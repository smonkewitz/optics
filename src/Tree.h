#pragma once

#include <sys/types.h>

#include <array>
#include <cstddef>
#include <limits>
#include <memory>

#include "Vec3.h"

namespace optics {

// A pointer-less node in a 3-d tree. A dimension, splitting value along that dimension,
// and the index of the point following the last point in the leaf is stored. The index
// of the first point in the node is obtained from the node to the left at the same
// level of the tree. Memory usage per node is just 16 bytes.
struct alignas(16) Node {
    static constexpr size_t SHIFT = 2;
    static constexpr size_t MASK = (static_cast<size_t>(1) << SHIFT) - 1;

    // Splitting value
    double split = std::numeric_limits<double>::quiet_NaN();
    // 2 LSBs: dimension of the splitting value (0, 1, 2), or 3 if the node is a leaf
    //   MSBs: index of first entry to the right of the split
    size_t metadata = static_cast<size_t>(-1);

    size_t right() const { return metadata >> SHIFT; }
    size_t splitDim() const { return metadata & MASK; }
    bool isLeaf() const { return splitDim() == 3; }

    void setSplitDim(size_t dim) { metadata = (dim & MASK) | (metadata & ~MASK); }
    void setRight(size_t index) { metadata = (index << SHIFT) | (metadata & MASK); }
};

constexpr size_t NOT_FOUND = static_cast<size_t>(-1);
constexpr size_t UNPROCESSED = static_cast<size_t>(-1);
constexpr size_t PROCESSED = static_cast<size_t>(-2);

// An entry in the data array to be indexed using a 3-d tree. It contains coordinates,
// along with the following:
//
// - An integer used to embed a singly linked list of range query results in the data
//   array.
// - A double used to store the distance of the point to the range query input point.
// - The reachability-distance of the point (defined by the OPTICS algorithm).
// - The CSV record from which coordinates were obtained.
//
// Memory usage per point is 64 bytes (a single cache line on most CPUs).
struct alignas(64) Point {
    // unit vector extracted from record
    Vec3 v;
    // distance to query point
    double dist = std::numeric_limits<double>::quiet_NaN();
    // OPTICS reachability distance
    double reach = std::numeric_limits<double>::infinity();
    // originating CSV record
    char const* record = nullptr;
    // index of next range query result or 0xffffffff
    size_t next = NOT_FOUND;
    // [UN]PROCESSED, or index in seed list
    size_t state = UNPROCESSED;
};

// A pointer-less 3-d tree class over an array of Point objects. Points belonging to a
// node are contiguous in memory. Furthermore, the location of the nodes themselves is
// implicit: the children of node i are located at positions 2*i + 1 and 2*i + 2 in an
// underlying array. Nodes therefore need not store pointers to their children, and
// siblings are contiguous in memory.
//
// The class supports a simple range query - finding all points within some distance D
// of a point. The result of this operation is returned as an index to the first Point
// in range - remaining results are available by traversal of the linked list embedded
// in the points. Because the results are expected to span a small number of k-d tree
// leaves and will already have been touched by the range query, the linked list is
// likely to be cache-resident prior to traversal. However, the consequence of this
// approach is that a tree and its associated Point array must only be used by a single
// thread at a time.
//
// It is also important to note that this class does not own the array of points over
// which it is defined - it is the caller's responsibility to ensure that the lifetime
// of the array exceeds the lifetime of the tree and that the array is not modified
// while the tree is alive.
class Tree {
   public:
    static constexpr size_t MAX_HEIGHT = sizeof(size_t) * 8 - 2;

    // Creates a new 3-d tree over an array of points. The tree construction process
    // modifies the order of points in the array but not the points themselves.
    //
    // - pointsPerLeaf:        Target # of points per leaf node, used to determine 3-d
    //                         tree height.
    // - leafExtentThreshold:  If the maximum extent of a 3-d tree node along each
    //                         dimension is below this number, then no children are
    //                         created for the node.
    Tree(Point* points, size_t numPoints, size_t pointsPerLeaf,
         double leafExtentThreshold);

    size_t size() const { return numPoints_; }
    size_t height() const { return height_; }
    Point const* getPoints() const { return points_; }

    // Locates all points in the 3-d tree within squared euclidian distance `dist` of
    // the input query point `v`.
    //
    // The result of the query is returned as a single integer index to the first point
    // in range - remaining results are available by traversal of the linked list
    // embedded in the point array. If no points are in range, NOT_FOUND is returned.
    size_t inRange(Vec3 const& v, double dist);

   private:
    Point* points_;  // unowned
    size_t numPoints_;
    size_t height_;
    std::unique_ptr<Node[]> nodes_;

    void build(double leafExtentThreshold);
};

}  // namespace optics
