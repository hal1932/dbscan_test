#include "dbscan.hpp"

#include <cstddef>
#include "vendor/octree/octree.h"

#include <type_traits>
#include <vector>
#include <iterator>
#include <memory>

// And this is the "dataset to kd-tree" adaptor class:

inline auto get_pt(const point2& p, std::size_t dim)
{
    if(dim == 0) return p.x;
    return p.y;
}


inline auto get_pt(const point3& p, std::size_t dim)
{
    if(dim == 0) return p.x;
    if(dim == 1) return p.y;
    return p.z;
}

inline auto get_pt(const sphere& p, std::size_t dim)
{
    if(dim == 0) return p.x;
    if(dim == 1) return p.y;
    if(dim == 2) return p.z;
    return p.radius;
}


template<typename Point>
struct PointAdapter
{
    const std::span<const Point>&  points;
    PointAdapter(const std::span<const Point>&  points) : points(points) { }

    /// CRTP helper method
    //inline const Derived& derived() const { return obj; }

    // Must return the number of data points
    inline std::size_t kdtree_get_point_count() const { return points.size(); }

    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline float kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
    {
        return get_pt(points[idx], dim);
    }

    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

    auto const * elem_ptr(const std::size_t idx) const
    {
        return &points[idx].x;
    }
};



auto sort_clusters(std::vector<std::vector<size_t>>& clusters)
{
    for(auto& cluster: clusters)
    {
        std::sort(cluster.begin(), cluster.end());
    }
}


using namespace OrthoTree;

auto search_range(const QuadtreeBoxC& space, const BoundingBox2D& box, float eps) -> std::vector<uint32_t>
{
    BoundingBox2D range = {
        box.Min[0] - eps, box.Min[1] - eps,
        box.Max[0] + eps, box.Max[1] + eps
    };
    return space.RangeSearch(range);
}

auto dbscan(const std::span<const aabb2>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>
{
    auto boxes = std::make_unique<BoundingBox2D[]>(data.size());
    for (size_t i = 0; i < data.size(); i++)
    {
        boxes[i] = {
            data[i].x_min, data[i].y_min,
            data[i].x_max, data[i].y_max,
        };
    }

    auto space = QuadtreeBoxC(std::span<BoundingBox2D>(boxes.get(), data.size()), 10);

    const auto n_points = data.size();
    auto visited = std::vector<bool>(n_points);
    auto clusters = std::vector<std::vector<size_t>>();

    for (size_t i = 0; i < n_points; i++)
    {
        if (visited[i])
        {
            continue;
        }

        auto matches = search_range(space, boxes[i], eps);
        if (matches.size() < static_cast<size_t>(min_pts))
        {
            continue;
        }
        visited[i] = true;

        auto cluster = std::vector({ i });

        while (!matches.empty())
        {
            auto nb_idx = matches.back();
            matches.pop_back();

            if (visited[nb_idx])
            {
                continue;
            }
            visited[nb_idx] = true;

            auto sub_matches = search_range(space, boxes[nb_idx], eps);
            if (sub_matches.size() >= static_cast<size_t>(min_pts))
            {
                std::copy(sub_matches.begin(), sub_matches.end(), std::back_inserter(matches));
            }

            //index.radiusSearch(adapt.elem_ptr(nb_idx), eps, sub_matches, SearchParams(32, 0.f, false));
            //if (sub_matches.size() >= static_cast<size_t>(min_pts))
            //{
            //    std::copy(sub_matches.begin(), sub_matches.end(), std::back_inserter(matches));
            //}

            cluster.push_back(nb_idx);
        }

        clusters.emplace_back(std::move(cluster));
    }

    sort_clusters(clusters);

    return clusters;
}
