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

            child_triangles.setConstant(-1);
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
         * if the triangle has children, child_triangles[i] is the index to the child triangle that
         * corconstructed with edge i
         */
        Eigen::Vector3i child_triangles;
        VectorT center;
        int parent = -1;
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
            set_as_parent(bt);
            auto polygon = get_surrouding_polygon(bt);

            reconstruct(polygon, i);
            update_parent_radius(bt);
        }

        // remove any triangle that has child or has vertex of the super triangle
        {
            int n = m_num_vertices;
            m_triangles.erase(
                std::remove_if(
                    m_triangles.begin(),
                    m_triangles.end(),
                    [n](const Triangle &t) {
                        return t.has_child() || t.is_super(n);
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

    Eigen::Matrix<double, 2, Eigen::Dynamic> get_voronoi_centers() const
    {
        int n = m_triangles.size();
        Eigen::Matrix<double, 2, Eigen::Dynamic> ret(2, n);
        for (int i = 0; i < n; ++i) {
            ret.col(i) = m_triangles[i].center;
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
        std::vector<std::set<int>> edges(m_num_vertices);

        auto insert = [&](const Triangle & t, int a, int b) {
            int p1 = t.point_indices(a);
            int p2 = t.point_indices(b);

            if (p1 > p2) {
                std::swap(p1, p2);
            }
            (void) edges[p1].insert(p2);
        };

        for (const auto & t : m_triangles) {
            insert(t, 0, 1);
            insert(t, 1, 2);
            insert(t, 2, 0);
        }

        // according to Euler's equation V-E+F = 1 for 2D plane
        int num_edges = m_num_vertices + m_triangles.size() -1;

        // convert set to matrix
        Eigen::Matrix<int, 2, Eigen::Dynamic> ret(2, num_edges);

        int i = 0;
        for (int p1 = 0; p1 < m_num_vertices; ++p1) {
            const auto & pnt_set = edges[p1];
            for (int p2 : pnt_set) {
                ret.col(i) << p1, p2;
                ++i;
            }
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
     * find all bad leaf triangles whose circumcircle covers the given point
     */
    std::set<int> find_bad_triangles(int point_index) const {
        std::set<int> triangles;
        find_bad_triangle_vertical_recur(0, point_index, triangles);
        return  triangles;
    }

    /*
     * search for a bad triangle that covers the point in the given triangle and its descendents
     */
    void find_bad_triangle_vertical_recur(int triangle_index, int point_index, std::set<int>& triangles) const
    {
        const Triangle & t = m_triangles[triangle_index];
        if (!t.circumcircle_covers(m_points.col(point_index))) {
            return;
        }

        if (t.has_child()) {
            for (int i = 0; i < 3; ++i) {
                int cti = t.child_triangles(i);
                if (cti > 0) {
                    find_bad_triangle_vertical_recur(cti, point_index, triangles);
                }
            }
        }
        else {
            triangles.insert(triangle_index);
        }
    }

    /*
     * get surrouding polygon based on collected bad triangles
     *
     * The result is matrix of indices, each column of it holds:
     * vertex and parent triangle (both are indices)
     *
     * Vertex[i] and vertex[i+1] forms edges i and the parent triangle is the triangle that formally
     * holds the edge. Note that edges belong to two bad triangles will be removed, so that for each
     * edge, there will be only one parent.
     */
    std::map<std::pair<int, int>, int> get_surrouding_polygon(const std::set<int> & triangle_index) const
    {
        auto reverse = [](const std::pair<int, int>& p) {
            return std::make_pair(p.second, p.first);
        };

        // map from edge to parent
        std::map<std::pair<int, int>, int> edge_records;

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
                    edge_records[edge] = ti;
                }
            }
        }

        return edge_records;
    }

    void set_as_parent(const std::set<int> & triangle_index)
    {
        for (int ti : triangle_index) {
            m_triangles[ti].is_parent = true;
        }
    }

    void reconstruct(const std::map<std::pair<int, int>, int> & polygon, int point_index)
    {
        for (const auto & edge : polygon) {
            int p1 = edge.first.first;
            int p2 = edge.first.second;
            int ti_parent = edge.second;

            int ti = (int)m_triangles.size();   // index to the new triangle

            m_triangles.emplace_back(m_points, p1, p2, point_index);
            Triangle & current_triangle = m_triangles.back();

            // set link from parent to this triangle
            current_triangle.parent = ti_parent;
            Triangle & parent_triangle = m_triangles[ti_parent];
            for (int j = 0; j < 3; ++j) {
                if (parent_triangle.point_indices[j] == p1) {
                    parent_triangle.child_triangles[j] = ti;
                    break;
                }
            }
        }
    }

    /*
     * update radius of the prarents so that they will cover all children
     */
    void update_parent_radius(const std::set<int> & ti_parent) {
        std::set<int> ti_grand_parent;
        for (int ti : ti_parent) {
            if (ti >= 0) {
                Triangle &t = m_triangles[ti];
                assert(t.is_parent);
                bool updated = false;
                for (int i= 0; i < 3; ++i) {
                    int ti_child = t.child_triangles(i);
                    const Triangle &t_child = m_triangles[ti_child];
                    if (ti_child >= 0) {
                        double dis_to_child = (t.center - t_child.center).norm();
                        if (t.radius < dis_to_child + t_child.radius) {
                            t.radius = dis_to_child + t_child.radius;
                            updated = true;
                        }
                    }
                }

                if (updated) {
                    ti_grand_parent.insert(t.parent);
                }
            }
        }

        if (!ti_grand_parent.empty()) {
            update_parent_radius(ti_grand_parent);
        }
    }

private:
    std::vector<Triangle> m_triangles;

    VerticesType m_points;
    int m_num_vertices;
};

