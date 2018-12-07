#pragma once

#include <algorithm>
#include <vector>
#include <map>
#include <set>

#include <Eigen/Geometry>

template<typename T>
class Delaunay
{
public:
    using VectorT = Eigen::Matrix<T, 2, 1>;
    using VerticesType = Eigen::Matrix<T, 2, Eigen::Dynamic>;

private:

    struct Triangle
    {
        Triangle(const VerticesType & points, int p1, int p2, int p3) :
            point_indices(p1, p2, p3)
        {
            Eigen::Matrix<T, 2, 2> D;
            D.col(0) = points.col(p2) - points.col(p1);
            D.col(1) = points.col(p3) - points.col(p1);

            VectorT U = 0.5 * D.transpose().inverse() *
                        D.colwise().squaredNorm().transpose();

            radius = U.norm();
            center = U + points.col(p1);

            linked_triangles.setConstant(-1);
        }

        bool circumcircle_covers(const VectorT & p) const
        {
            return (p - center).norm() <= radius;
        }

        bool has_child() const
        {
            return is_parent;
        }

        bool is_super(int n) const
        {
            return (point_indices.array() >= n).any();
        }

        Eigen::Vector3i point_indices;

        /*
         * Let edge i of the triangle be:
         * edge 0: (point_indices[0], point_indices[1])
         * edge 1: (point_indices[1], point_indices[2])
         * edge 2: (point_indices[2], point_indices[0])
         *
         * linked triangles are
         * 1. child trianlges, if the triangle is a parent triangle
         *    in this case, linked_triangles[i] is constructed with edge i
         * 2. neighbor triangles, if the triangle is a leaf triangle
         *    in this case, linked_triangles[i] shares the edge i
         */
        Eigen::Vector3i linked_triangles;
        VectorT center;
        double radius;
        bool is_parent = false;
    };

public:
    /**
     * Constructor
     * @param vertices
     * The points to be trianglized. Each column is a point
     */
    Delaunay(const VerticesType & points) : m_points(2, points.cols()+3),
            m_num_vertices(points.cols())
    {
        m_points.leftCols(m_num_vertices) = points;
    }

    /**
     * Start triangulation
     */
    void triangulate()
    {
        create_super_triangle();

        // trianglate with new points
        for (int i = 0; i < m_num_vertices; ++i) {

            std::set<int> bt = find_bad_triangles(i);
            Eigen::MatrixXi polygon = get_surrouding_polygon(bt);
            set_as_parent(bt);

            reconstruct(polygon, i);
        }

        // remove any triangle that has child or has vertex of the super triangle
        {
            int n = m_num_vertices;
            m_triangles.erase(
                std::remove_if(
                    m_triangles.begin(),
                    m_triangles.end(),
                    [n](const Triangle &t) {
                        return t.is_super(n) || t.has_child();
                    }
                ),
                m_triangles.end()
            );
        }
    }

    /**
     * Get triangles after triangulation
     *
     * @return
     * each column holds indices to the three points in the original input
     */
    Eigen::Matrix<int, 3, Eigen::Dynamic> get_triangles() const
    {
        int n = m_triangles.size();
        Eigen::Matrix<int, 3, Eigen::Dynamic> ret(3, n);
        for (int i = 0; i < n; ++i) {
            ret.col(i) = m_triangles[i].point_indices;
        }
        return ret;
    }

    /**
     * Get triangles after triangulation
     *
     * @return
     * each column holds indices to the two points in the original input
     */
    Eigen::Matrix<int, 2, Eigen::Dynamic> get_edges() const
    {
        std::set<std::pair<int, int>> edges;
        auto insert = [&](const Triangle & t, int a, int b) {
            int p1 = t.point_indices(a);
            int p2 = t.point_indices(b);

            if (p1 > p2) {
                std::swap(p1, p2);
            }
            (void) edges.insert(std::make_pair(p1, p2));
        };

        for (const auto & t : m_triangles) {
            insert(t, 0, 1);
            insert(t, 1, 2);
            insert(t, 2, 0);
        }

        // convert set to matrix
        int n = edges.size();
        Eigen::Matrix<int, 2, Eigen::Dynamic> ret(2, n);
        int i = 0;
        for (const auto & item : edges) {
            ret.col(i) << item.first, item.second;
            ++i;
        }

        return ret;
    }

private:

    void create_super_triangle()
    {
        VectorT minp = m_points.leftCols(m_num_vertices).rowwise().minCoeff();
        VectorT maxp = m_points.leftCols(m_num_vertices).rowwise().maxCoeff();
        VectorT dp = maxp - minp;
        VectorT center = (minp + maxp) / 2;

        dp *= 10;    // make it much larger than the scale of the point set
                     // it is crucial to make the super triangle removable

        // constructing the super triangle vertices into point list

        Eigen::Rotation2D<T> rotation(2.0/3.0 * EIGEN_PI);

        int n = m_num_vertices;
        m_points.col(n) = center + dp;
        m_points.col(n+1) = center + rotation * dp;
        m_points.col(n+2) = center + rotation.inverse() * dp;

        // add the supper triangle into triangle list
        m_triangles.emplace_back(m_points, n, n+1, n+2);
    }

    /*
     * find all bad triangles whose circumcircle covers the given point
     *
     * It searches the bad triangle tree rooted by the super triangle to get one bad leaf triangle
     * and then finds the leaf bad triangles by searching the neighbors of the first-found bad
     * leaf triangle
     */
    std::set<int> find_bad_triangles(int point_index) const {
        std::set<int> triangles;

        int bt;
        {
            // find a bad triangle vetically
            auto result = find_bad_triangle_recur(0, point_index);
            assert(result.first);

            bt = result.second;
        }

        // then find bad triangles horizontally
        find_bad_triangles(bt, point_index, triangles);

        return  triangles;
    }

    /*
     * search for a bad triangle that covers the point in the given triangle and its descendents
     */
    std::pair<bool, int> find_bad_triangle_recur(int triangle_index, int point_index) const
    {
        const Triangle & t = m_triangles[triangle_index];

        if (t.has_child()) {
            for (int i = 0; i < 3; ++i) {
                int cti = t.linked_triangles(i);
                if (cti > 0) {
                    auto result = find_bad_triangle_recur(cti, point_index);
                    if (result.first) {
                        return  result;
                    }
                }
            }
        }
        else {
            if (t.circumcircle_covers(m_points.col(point_index))) {
                return std::make_pair(true, triangle_index);
            }
        }

        return std::make_pair(false, -1);
    }

    // find leaf bad triangles by searching neighbors
    void find_bad_triangles(int triangle_index, int point_index, std::set<int>& triangles) const
    {
        if (triangle_index == -1) {
            return;
        }

        auto it = triangles.find(triangle_index);
        if (it != triangles.end()) {
            // has been added, ignore it
            return;
        }

        const Triangle & t = m_triangles[triangle_index];
        if (!t.circumcircle_covers(m_points.col(point_index))) {
            return;
        }

        triangles.insert(triangle_index);

        for (int i = 0; i < 3; ++i) {
            find_bad_triangles(t.linked_triangles(i), point_index, triangles);
        }
    }

    /*
     * get surrouding polygon based on collected bad triangles
     *
     * The result is matrix of indices, each column of it holds:
     * vertex, parent triangle and neighbor triangle
     *
     * Vertex[i] and vertex[i+1] forms edges i and the parent triangle is the triangle that formally
     * holds the edge. Note that edges belong to two bad triangles will be removed, so that for each
     * edge, there will be only one parent.
     *
     * The neighor triangle for the edge is the 'normal' triangle that shares the edge. This is to
     * be used to fill in the neighbor record for the new triangles.
     */
    Eigen::MatrixXi get_surrouding_polygon(const std::set<int> & triangle_index) const
    {
        auto reverse = [](const std::pair<int, int>& p) {
            return std::make_pair(p.second, p.first);
        };

        // from edge to (parent, neighbor)
        std::map<std::pair<int, int>, std::pair<int, int>> edge_records;

        for (int ti : triangle_index) {
            const Triangle & t = m_triangles[ti];

            for (int i = 0; i < 3; ++i) {
                int j = (i+1) % 3;
                auto edge = std::make_pair(t.point_indices(i), t.point_indices(j));

                auto it = edge_records.find(reverse(edge));
                if (it != edge_records.end()) {
                    // this edge was added by another bad triangle, that means this edge is internal
                    // it cannot be used to contruct new triangles.
                    edge_records.erase(it);
                }
                else {
                    int neighbor = t.linked_triangles(i);
                    edge_records[edge] = std::make_pair(ti, neighbor);
                }
            }
        }

        // convert the map to matrix
        int n = (int)edge_records.size();
        Eigen::MatrixXi ret(3, n);
        auto it = edge_records.begin();
        int i = 0;
        do {
            ret(0, i) = it->first.first;
            ret(1, i) = it->second.first;
            ret(2, i) = it->second.second;

            it = edge_records.lower_bound(std::make_pair(it->first.second, 0));
            ++i;
        } while (it != edge_records.begin() && it != edge_records.end());

        assert(i == n);

        return ret;
    }

    void set_as_parent(const std::set<int> & triangle_index) {
        for (int ti : triangle_index) {
            Triangle &t = m_triangles[ti];
            t.is_parent = true;
            t.linked_triangles.setConstant(-1);
        }
    }

    void reconstruct(const Eigen::MatrixXi & polygon, int point_index)
    {
        int ti_base = (int)m_triangles.size();  // index to the first added triangle

        int n = polygon.cols();
        for (int i = 0; i < n; ++i) {
            int p1 = polygon(0, i);
            int j = (i+1) % n;
            int p2 = polygon(0, j);

            int ti = (int)m_triangles.size();   // index to the new triangle

            m_triangles.emplace_back(m_points, p1, p2, point_index);
            Triangle & current_triangle = m_triangles.back();

            // set link from parent to this triangle
            Triangle & parent_triangle = m_triangles[polygon(1, i)];
            for (int j = 0; j < 3; ++j) {
                if (parent_triangle.point_indices[j] == p1) {
                    parent_triangle.linked_triangles[j] = ti;
                    break;
                }
            }

            // set neighbor
            int ti_neighbor = polygon(2, i);
            if (ti_neighbor > 0) {
                Triangle &neighbor_triangle = m_triangles[ti_neighbor];
                for (int j = 0; j < 3; ++j) {
                    // the neighbor triangle has the edge p2->p1
                    if (neighbor_triangle.point_indices[j] == p2) {
                        neighbor_triangle.linked_triangles[j] = ti;

                        break;
                    }
                }
            }

            current_triangle.linked_triangles <<
                    ti_neighbor,
                    ((i == (n-1))? ti_base :ti + 1),
                    ((i == 0)? (ti_base + n - 1): (ti - 1));
        }
    }

private:
    std::vector<Triangle> m_triangles;

    VerticesType m_points;
    int m_num_vertices;
};

