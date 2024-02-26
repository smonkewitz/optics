#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "Tree.h"
#include "Vec3.h"

namespace optics {
namespace {

struct MatchOracle {
    Vec3 query;
    std::vector<size_t> expectedMatches;

    bool isExpected(size_t match) const {
        return std::find(expectedMatches.cbegin(), expectedMatches.cend(), match) !=
               expectedMatches.cend();
    }
};

constexpr double TestRadius = 0.25;

void MakeTestPoints(std::vector<Point>& points, std::vector<MatchOracle>& queries) {
    std::mt19937_64 rng(1234);
    double const deltaDec = 4.0 * TestRadius;

    size_t pointId = 0;
    // divide the unit sphere into longitude/latitude angle boxes
    for (double dec = -90.0; dec < 90.0; dec += deltaDec) {
        double const maxAbsDec = std::max(std::fabs(dec), std::fabs(dec + deltaDec));
        double const raExtent = LongitudeExtent(4.0 * TestRadius, maxAbsDec);
        size_t const numBoxes = static_cast<size_t>(std::floor(360.0 / raExtent));
        double const deltaRa = 360.0 / numBoxes;
        for (size_t i = 0; i < numBoxes; ++i) {
            double const ra = i * deltaRa;
            // create a random query point inside a sub region of each box such that a
            // circle of the given radius centered on that point is guaranteed not to
            // cross the box boundaries
            MatchOracle oracle;
            auto queryLonLat =
                LonLat::random(rng, ra + deltaRa * 0.4, ra + deltaRa * 0.6,
                               dec + deltaDec * 0.4, dec + deltaDec * 0.6);
            oracle.query = queryLonLat;
            int const numPointsToGenerate =
                std::uniform_int_distribution<int>{0, 64}(rng);
            for (int j = 0; j < numPointsToGenerate; ++j) {
                auto perturbedLonLat = queryLonLat.perturb(rng, TestRadius);
                double distance = queryLonLat.distance(perturbedLonLat);
                if (distance >= 1.5 * TestRadius) {
                    continue;
                }
                Point point;
                point.v = perturbedLonLat;
                if (distance < 0.999999 * TestRadius) {
                    // perturbedLonLat is definitely a match for queryLonLat
                    point.state = pointId;
                    oracle.expectedMatches.push_back(pointId);
                    points.push_back(point);
                    ++pointId;
                } else if (distance > 1.0000001 * TestRadius) {
                    // perturbedLonLat is definitely not a match for queryLonLat
                    point.state = pointId;
                    points.push_back(point);
                    ++pointId;
                }
            }
            queries.push_back(oracle);
        }
    }
}

TEST(TreeTest, InRange) {
    std::vector<Point> points;
    std::vector<MatchOracle> queries;

    double const distance = SquaredEuclidianDistance(TestRadius);
    MakeTestPoints(points, queries);
    Tree tree{points.data(), points.size(), 32, 0.0};
    std::vector<size_t> matches;
    for (auto const& oracle : queries) {
        size_t index = tree.inRange(oracle.query, distance);
        matches.clear();
        while (index != static_cast<size_t>(-1)) {
            matches.push_back(points[index].state);
            index = points[index].next;
        }
        EXPECT_THAT(matches,
                    testing::UnorderedElementsAreArray(oracle.expectedMatches));
    }
}

}  // namespace
}  // namespace optics
