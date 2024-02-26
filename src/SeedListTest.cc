#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <cstddef>
#include <limits>
#include <random>
#include <vector>

#include "SeedList.h"
#include "Tree.h"

namespace optics {
namespace {

std::vector<Point> const MakePoints(size_t n) {
    std::vector<Point> points{n, Point{}};
    for (size_t i = 0; i < n; ++i) {
        points[i].reach = i;
    }
    return points;
}

std::vector<Point> const MakePoints(size_t n, std::mt19937_64& rng) {
    std::vector<Point> points{n, Point{}};
    for (size_t i = 0; i < n; ++i) {
        points[i].reach = std::uniform_int_distribution<size_t>{0, n >> 1}(rng);
    }
    return points;
}

// Tests add() and pop() methods of SeedList class
TEST(SeedListTest, AddPopBasic) {
    size_t n = 128;
    // construct points with strictly increasing reachability distance
    auto points = MakePoints(n);
    SeedList sl(points.data(), n);
    EXPECT_TRUE(sl.empty());
    EXPECT_EQ(sl.capacity(), n);
    EXPECT_EQ(sl.pop(), static_cast<size_t>(-1));
    EXPECT_TRUE(sl.checkInvariants());
    sl.add(0);
    EXPECT_TRUE(sl.checkInvariants());
    EXPECT_EQ(sl.pop(), 0);
    EXPECT_EQ(sl.size(), 0);
    EXPECT_TRUE(sl.checkInvariants());
    sl.add(n - 1);
    sl.add(0);
    EXPECT_TRUE(sl.checkInvariants());
    EXPECT_EQ(sl.pop(), 0);
    EXPECT_TRUE(sl.checkInvariants());
    EXPECT_EQ(sl.pop(), n - 1);
    EXPECT_EQ(sl.size(), 0);
    // add points in increasing reachability-distance order
    for (size_t i = 0; i < n; ++i) {
        sl.add(i);
    }
    EXPECT_EQ(sl.size(), n);
    EXPECT_TRUE(sl.checkInvariants());
    // check that points are popped in increasing reachability-distance order
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(sl.pop(), i);
        EXPECT_TRUE(sl.checkInvariants());
    }
    EXPECT_EQ(sl.size(), 0);
    // add points in decreasing reachability-distance order
    for (size_t i = n; i > 0; --i) {
        sl.add(i - 1);
    }
    EXPECT_EQ(sl.size(), n);
    EXPECT_TRUE(sl.checkInvariants());
    // check that points are popped in increasing reachability-distance order
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(sl.pop(), i);
        EXPECT_TRUE(sl.checkInvariants());
    }
}

// Tests add() and pop() methods of SeedList with randomly ordered inputs
TEST(SeedListTest, AddPopRandom) {
    std::mt19937_64 rng(1234);
    size_t n = 127;
    auto points = MakePoints(n, rng);
    SeedList sl(points.data(), n);

    for (size_t i = 0; i < n; ++i) {
        sl.add(i);
        EXPECT_TRUE(sl.checkInvariants());
    }
    EXPECT_EQ(sl.size(), n);
    double maxReach = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i) {
        double reach = points[sl.pop()].reach;
        EXPECT_TRUE(sl.checkInvariants());
        EXPECT_TRUE(reach >= maxReach);
        maxReach = reach;
    }
    EXPECT_EQ(sl.size(), 0);
}

// Tests the update() method of SeedList
TEST(SeedListTest, Update) {
    size_t n = 120;
    auto points = MakePoints(n);
    std::vector<size_t> order;
    SeedList sl(points.data(), n);
    for (size_t i = 0; i < n; ++i) {
        sl.add(i);
    }
    for (size_t i = 0; i < n; ++i) {
        order.push_back(sl.pop());
    }
    for (size_t i = 0; i < n; ++i) {
        sl.add(i);
    }
    for (size_t i = 0; i < n; ++i) {
        sl.update(i, -points[i].reach);
        EXPECT_TRUE(sl.checkInvariants());
    }
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(sl.pop(), order[n - i - 1]);
    }
}

}  // namespace
}  // namespace optics
