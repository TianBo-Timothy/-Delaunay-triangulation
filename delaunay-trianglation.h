#pragma once

#include <vector>
#include <algorithm>

#include <Eigen/Dense>

template<typename T>
class Delaunay
{
public:
    using VectorT = Eigen::Matrix<T, 2, 1>;
    using VerticesType = Eigen::Matrix<T, 2, Eigen::Dynamic>;

private:
    struct Edge
    {
        Edge(int _p1, int _p2) :
            p1(_p1), p2(_p2)
        {
            if (p1 > p2) {
                std::swap(p1, p2);
            }
        }

        bool operator<(const Edge & rhs) {
            return std::make_pair(p1, p2) < std::make_pair(rhs.p1, rhs.p2);
        }

        int p1;
        int p2;
    };

    struct Triangle
    {
        Triangle(const VerticesType & vertices, int p1, int p2, int p3) :
            point_indices(p1, p2, p3)
        {
            Matrix<T, 2, 2> D;
            D.col(0) = vertices.col(p2) - vertices.col(p1);
            D.col(1) = vertices.col(p3) - vertices.col(p1);

            VectorT U = D.inverse() * D.colwise().squaredNorm().transpose();
            radius = U.norm();
            center = U + vertices.col(p1);
        }

        bool has(int point_index) const
        {
            return (point_indices.array() == point_index).any();
        }

        bool circumcircle_covers(const VectorT & p) const
        {
            return (p - center).norm() <= radius;
        }

        void invalidate()
        {
            point_indices(0) = -1;
        }

        void is_bad() const
        {
            return (point_indices(0) == -1);
        }


        Vector3i point_indices;
        VectorT center;
        double radius;
    };

public:
    /**
     * Constructor
     * @param vertices
     * The points to be trianglized. Each column is a point
     */
    Delaunay(const VerticesType & vertices) : m_vertices(vertices),
            m_num_vertices(vertices.cols())
    {
    }

    /**
     * Start triangulation
     */
    void triangulate()
    {
        VectorT minp = m_vertices.rowwise().minCoef();
        VectorT maxp = m_vertices.rowwise().maxCoef();
        VectorT dp = maxp - minp;
        VectorT midp = (minp + maxp) / 2;

        // constructing the super triangle
        VectorT p1, p2, p3;
        {
            Matrix<T, 2, 2> rotate_mat;
            const double theta = 3.1415926535 / 3;
            rotate_mat << std::cos(theta), std::sin(theta),
                         -std::sin(theta), std::cos(theta);
            p1 = midp + dp;
            p2 = midp + rotate_mat * dp;
            p3 = midp + rotate_mat.transpose() * dp;
        }

        // add vertices of the super triangle into vertices
        {
            int n = m_num_vertices;
            m_vertices.conservetiveResize(2, n+3);
            m_vertices.col(n) = p1;
            m_vertices.col(n+1) = p2;
            m_vertices.col(n+2) = p3;
        }

        // add the supper triangle into triangle list
        {
            int n = m_num_vertices;
            m_triangles.push_back({n, n+1, n+2});
        }

        // trianglate with new points
        for (int i = 0; i < m_num_vertices; ++i) {

            std::set<Edge> polygon;
            for (auto & t : m_triangles) {
                if (t.circumcircle_covers(m_vertices.col(i))) {
                    add_triangle_into_polygon(polygon, t);
                    t.invalidate();
                }
            }

            // remove 'bad' triangles
            m_triangles.erase(
                std::remove_if(
                    begin(m_triangles),
                    end(m_triangles),
                    [](const Triangle &t){ return t.is_bad(); }
                ),
                end(m_triangles)
            );

            // create new triangles with point i and the polygon
            for(const auto & e : polygon) {
                m_triangles.push_back(Triangle(e.p1, e.p2, i));
            }
        }

        // remove any triangle that has vertex of the super triangle
        m_triangles.erase(
            std::remove_if(
                m_triangles.begin(),
                m_triangles.end(),
                [](const Triangle &t) {
                    return t.has(m_num_vertices) || t.has(m_num_vertices+1) ||
                            t.has(m_num_vertices+2);
                }
            ),
            m_triangles.end()
        );

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
            ret.col(i) = m_triangles.point_indices;
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
        std::set<Edge> & edges;
        for (auto & t : m_triangles) {
            (void)edges.insert(Edge(t.point_indices(0), t.point_indices(1)));
            (void)edges.insert(Edge(t.point_indices(1), t.point_indices(2)));
            (void)edges.insert(Edge(t.point_indices(2), t.point_indices(0)));
        }

        int n = edges.size();
        Eigen::Matrix<int, 2, Eigen::Dynamic> ret(2, n);
        int i = 0;
        for (const auto & item : edges) {
            ret.col(i) << item.p1, item.p2;
        }

        return ret;
    }

private:
    // add edges of the triangle into the polygon
    void add_triangle_into_polygon(std::set<Edge> & ploygon, const Triangle & t)
    {
        add_edge_into_polygon(Edge(t.point_indices(0), t.point_indices(1)));
        add_edge_into_polygon(Edge(t.point_indices(1), t.point_indices(2)));
        add_edge_into_polygon(Edge(t.point_indices(0), t.point_indices(2)));
    }

    // add edge into the polygon
    // remove the edge if it already exists
    void add_edge_into_polygon(std::set<Edge> & ploygon, const Edge & e)
    {
        auto result = polygon.insert(e);
        if (result.second == false) {
            polygon.erase(result.first);
        }
    }

private:
    std::vector<Triangle> m_triangles;

    VerticesType m_vertices;
    int m_num_vertices;
};

