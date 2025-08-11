#pragma once


#include <cassert>
#include <cstddef>
#include <span>
#include <vector>
#include <cstdlib>

struct point2
{
    float x, y;
};

struct point3
{
    float x, y, z;
};

auto dbscan(const std::span<const point2>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>;
auto dbscan(const std::span<const point3>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>;

struct sphere
{
    float x, y, z;
    float radius;
};

auto dbscan(const std::span<const sphere>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>;

struct aabb2
{
    float x_min, y_min;
    float x_max, y_max;
    float z;
};

struct aabb3
{
    float x_min, y_min, z_min;
    float x_max, y_max, z_max;
};

auto dbscan(const std::span<const aabb2>& data, float eps, int min_pts) -> std::vector<std::vector<size_t>>;

// template<size_t dim>
// auto dbscan(const std::span<float>& data, float eps, int min_pts)
// {
//     static_assert(dim == 2 or dim == 3, "This only supports either 2D or 3D points");
//     assert(data.size() % dim == 0);
    
//     if(dim == 2)
//     {
//         auto * const ptr = reinterpret_cast<float const*> (data.data());
//         auto points = std::span<const point2> 
//     }
// }
