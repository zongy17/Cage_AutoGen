#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__


#include <iostream>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <vector>

#include <Eigen/Dense>

typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector3d Point3d;

class Triangle 
{
public:
    Triangle(Point3d pt0, Point3d pt1, Point3d pt2) : m_pt {pt0, pt1, pt2} 
    {
        m_pt[0] = pt0; m_pt[1] = pt1; m_pt[2] = pt2;
        segs[0] = m_pt[1] - m_pt[0];
        segs[1] = m_pt[2] - m_pt[1];
        segs[2] = m_pt[0] - m_pt[2];
        m_vecNormal = segs[0].cross(segs[1]);
        m_vecNormal.normalize();// 单位化法向量
    }

    // 有正负之分
    double GetDistanceFromPointToTrianglePlane(Point3d pt) const {
        auto Pt0PtT = pt - m_pt[0];// 从三角形的点0指向pt
        return m_vecNormal.dot(Pt0PtT);
    }

    void GetNormal(Eigen::Vector3d &vecNormal) const {
        vecNormal = m_vecNormal;
    }

    Point3d m_pt[3];
    vec3d segs[3];
    Eigen::Vector3d m_vecNormal;
};


#define EPSION 1e-7

enum SegmentTriangleType {
    LYON,       // 刚好躺在上面(2个交点)
    HITTHROUGH, // 线段穿过三角形面片(1个交点)
    PASSBY      // 线段没有穿过(0个交点)
};

enum PointTriangleType {
    IN, // 点位于三角形内
    OUT // 点位于三角形外
};

enum TriangleTriangleType {
    INTERSECT, // 相交(2个交点)
    TOUCH,     // 刚好某个点挨着另一个三角形（1个交点）
    SEGREGATE  // 相离(0个交点)
};

bool IsZero(double value, double epsion = EPSION);

bool IsEqual(double v1, double v2, double epsion = EPSION);

bool IsPositive(double value, double epsion = EPSION);

bool IsNegative(double value, double epsion = EPSION);

int GetSignType(double value);

PointTriangleType testPointTriangle(Point3d const & pt, Triangle const & tri);

SegmentTriangleType testSegmentTriangle(Point3d const& from, Point3d const& to, Triangle const& tri, Point3d & hitPoint);

TriangleTriangleType testTriangleTriangle(Triangle const& tri1, Triangle const& tri2, std::vector<Point3d> & hitPoints);

#endif