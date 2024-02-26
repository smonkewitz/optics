#include "SeedList.h"

#include <absl/log/check.h>
#include <absl/log/log.h>

#include <cstdint>

#include "Tree.h"

namespace optics {

size_t SeedList::pop() {
    size_t s = size_;
    if (s == 0) {
        return NOT_FOUND;
    }
    size_t smallest = heap_[0];
    points_[smallest].state = PROCESSED;
    size_ = --s;
    if (s > 1) {
        siftDown(heap_[s]);
    } else if (s == 1) {
        size_t i = heap_[1];
        heap_[0] = i;
        points_[i].state = 0;
    }
    return smallest;
}

void SeedList::add(size_t i) {
    DCHECK(i < capacity());
    DCHECK(size() < capacity());
    size_t s = size_;
    size_ = s + 1;
    if (s == 0) {
        heap_[0] = i;
        points_[i].state = 0;
    } else {
        siftUp(s, i);
    }
}

void SeedList::update(size_t i, double reach) {
    DCHECK(i < capacity());
    size_t heapIndex = points_[i].state;
    if (heapIndex < PROCESSED) {
        DCHECK(heap_[heapIndex] == i);
        // the i-th point is already in the seed list
        if (reach < points_[i].reach) {
            points_[i].reach = reach;
            siftUp(heapIndex, i);
        }
    } else {
        // add i-th point to the seed list
        points_[i].reach = reach;
        add(i);
    }
}

void SeedList::siftUp(size_t heapIndex, size_t pointIndex) {
    DCHECK(heapIndex < size());
    DCHECK(pointIndex < capacity());
    double reach = points_[pointIndex].reach;
    while (heapIndex > 0) {
        size_t parentHeapIndex = (heapIndex - 1) >> 1;
        size_t parentPointIndex = heap_[parentHeapIndex];
        if (points_[parentPointIndex].reach <= reach) {
            break;
        }
        heap_[heapIndex] = parentPointIndex;
        points_[parentPointIndex].state = heapIndex;
        heapIndex = parentHeapIndex;
    }
    heap_[heapIndex] = pointIndex;
    points_[pointIndex].state = heapIndex;
}

void SeedList::siftDown(size_t pointIndex) {
    DCHECK(pointIndex < capacity());
    double reach = points_[pointIndex].reach;
    size_t halfSize = size_ >> 1;
    size_t heapIndex = 0;
    while (heapIndex < halfSize) {
        size_t childHeapIndex = (heapIndex << 1) + 1;
        size_t siblingHeapIndex = childHeapIndex + 1;
        size_t childPointIndex = heap_[childHeapIndex];
        double childReach = points_[childPointIndex].reach;
        if (siblingHeapIndex < size_) {
            size_t siblingPointIndex = heap_[siblingHeapIndex];
            double siblingReach = points_[siblingPointIndex].reach;
            if (siblingReach < childReach) {
                childReach = siblingReach;
                childPointIndex = siblingPointIndex;
                childHeapIndex = siblingHeapIndex;
            }
        }
        if (reach <= childReach) {
            break;
        }
        heap_[heapIndex] = childPointIndex;
        points_[childPointIndex].state = heapIndex;
        heapIndex = childHeapIndex;
    }
    heap_[heapIndex] = pointIndex;
    points_[pointIndex].state = heapIndex;
}

bool SeedList::checkInvariants() const {
    // check that each point knows its location in the seed list
    for (size_t i = 0; i < numPoints_; ++i) {
        size_t h = points_[i].state;
        if (h < PROCESSED) {
            if (h >= size_) {
                LOG(ERROR) << "point " << i << " has invalid seed list index " << h
                           << " >= " << size_;
                return false;
            }
            if (heap_[h] != i) {
                LOG(ERROR) << "point " << i << " has an incorrect seed list index "
                           << h;
                return false;
            }
        }
    }
    for (size_t i = 0; i < size_; ++i) {
        size_t p = heap_[i];
        if (p >= numPoints_) {
            LOG(ERROR) << "heap entry " << i << " has invalid point index " << p
                       << " >= " << numPoints_;
            return false;
        }
        if (points_[p].state != i) {
            LOG(ERROR) << "point " << p << " has incorrect seed list index != " << i;
            return false;
        }
    }
    // check the heap invariant
    for (size_t i = 0; i < size_ >> 1; ++i) {
        double reach = points_[heap_[i]].reach;
        size_t h = (i << 1) + 1;
        size_t p = heap_[h];
        if (points_[p].reach < reach) {
            LOG(ERROR) << "heap invariant violation";
            return false;
        }
        h += 1;
        if (h < size_) {
            size_t p = heap_[h];
            if (points_[p].reach < reach) {
                LOG(ERROR) << "heap invariant violation";
                return false;
            }
        }
    }
    return true;
}

}  // namespace optics
