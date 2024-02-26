#pragma once

#include <cstddef>
#include <memory>

#include "Tree.h"

namespace optics {

// Class for managing the OPTICS ordered seeds. Methods to add a seed, remove the seed
// with the smallest reachability distance and to decrease the reachability of a seed
// are provided.
class SeedList {
   public:
    SeedList(Point* points, size_t numPoints)
        : heap_{std::make_unique<size_t[]>(numPoints)},
          points_{points},
          size_{0},
          numPoints_{numPoints} {}

    bool empty() const { return size_ == 0; }
    size_t size() const { return size_; }
    size_t capacity() const { return numPoints_; }

    // Finds the point with the smallest reachability-distance, removes it from the seed
    // list, and returns its index. If the seed list is empty, returns NOT_FOUND.
    size_t pop();

    // Adds the i-th point to this seed list. Assumes that:
    // - i < capacity()
    // - size() < capacity()
    void add(size_t i);

    // Updates the reachability-distance of the i-th point. If it isn't already in the
    // seed list, it is added. Otherwise, if the new reachability-distance is smaller
    // than the current one, the i-th point's reachability-distance is updated. Assumes
    // that i < capacity().
    void update(size_t i, double reach);

    bool checkInvariants() const;

   private:
    std::unique_ptr<size_t[]> heap_;
    Point* points_;  // unowned
    size_t size_;
    size_t numPoints_;

    void siftUp(size_t heapIndex, size_t pointIndex);
    void siftDown(size_t pointIndex);
};

}  // namespace optics
