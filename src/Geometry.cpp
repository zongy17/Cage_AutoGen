#include "Geometry.h"

bool IsZero(double value, double epsion) {
    return std::abs(value) < epsion;
}

bool IsEqual(double v1, double v2, double epsion) {
    return IsZero(v1-v2, epsion);
}

bool IsPositive(double value, double epsion) {
    return value > epsion;
}

bool IsNegative(double value, double epsion) {
    return value < - epsion;
}

int GetSignType(double value)
{
    if (IsZero(value)) return 0;
    if (IsPositive(value)) return 1;
    return -1;
}

bool isZeorVec3(vec3d const & vec) {
    // return IsZero(vec.x()) && IsZero(vec.y()) && IsZero(vec.z());
    return IsZero(vec.norm());
}

PointTriangleType testPointTriangle(Point3d const & pt, Triangle const & tri) {
    vec3d vec0 = tri.segs[0].cross(pt - tri.m_pt[0]); if (isZeorVec3(vec0)) return IN;// 刚好落在三角形的边上
    vec3d vec1 = tri.segs[1].cross(pt - tri.m_pt[1]); if (isZeorVec3(vec1)) return IN;// 刚好落在三角形的边上
    vec3d vec2 = tri.segs[2].cross(pt - tri.m_pt[2]); if (isZeorVec3(vec2)) return IN;// 刚好落在三角形的边上
    vec0.normalize();
    vec1.normalize();
    vec2.normalize();
    double dir01 = vec0.dot(vec1);  int sign_dir01 = GetSignType(dir01);
    double dir12 = vec1.dot(vec2);  int sign_dir12 = GetSignType(dir12);
    double dir20 = vec2.dot(vec0);  int sign_dir20 = GetSignType(dir20);
    // assert(sign_dir01 != 0);
    // assert(sign_dir12 != 0);
    // assert(sign_dir20 != 0);
    if (sign_dir01==0 || sign_dir12==0 || sign_dir20==0) {
        printf(" pt: %.12lf %.12lf %.12lf\n", pt(0), pt(1), pt(2));
        printf(" tri 0: %.12lf %.12lf %.12lf\n", tri.m_pt[0](0), tri.m_pt[0](1), tri.m_pt[0](2));
        printf(" tri 1: %.12lf %.12lf %.12lf\n", tri.m_pt[1](0), tri.m_pt[1](1), tri.m_pt[1](2));
        printf(" tri 2: %.12lf %.12lf %.12lf\n", tri.m_pt[2](0), tri.m_pt[2](1), tri.m_pt[2](2));
        printf(" vec0: %.12lf %.12lf %.12lf\n", vec0(0), vec0(1), vec0(2));
        printf(" vec1: %.12lf %.12lf %.12lf\n", vec1(0), vec1(1), vec1(2));
        printf(" vec2: %.12lf %.12lf %.12lf\n", vec2(0), vec2(1), vec2(2));
        printf(" dir01: %.12lf, dir12: %.12lf, dir20: %.12lf\n", dir01, dir12, dir20);
        exit(1);
    }
    
    if (sign_dir01==sign_dir12 && sign_dir12==sign_dir20)// pt在三条边的同一侧
        return IN;
    else 
        return OUT;
}

SegmentTriangleType testSegmentTriangle(Point3d const& from, Point3d const& to, Triangle const& tri, Point3d & hitPoint) {
    vec3d seg = to - from;
    double dist_from = tri.GetDistanceFromPointToTrianglePlane(from); int sign_from = GetSignType(dist_from);
    double dist_to   = tri.GetDistanceFromPointToTrianglePlane(to);   int sign_to   = GetSignType(dist_to);
    if (sign_from == sign_to){
        if (sign_from == 0) 
            return LYON;// 2个交点
        else
            return PASSBY;// 0个交点
    } else if (sign_from == 0) {// from位于tri的平面上
        if (testPointTriangle(from, tri) == IN){// from刚好位于tri内
            hitPoint = from;
            return HITTHROUGH;
        } else // from在tri外，则不可能有交点
            return PASSBY;
    } else if (sign_to   == 0) {
        if (testPointTriangle(to, tri) == IN){// to刚好位于tri内
            hitPoint = to;
            return HITTHROUGH;
        } else // to在tri外，则不可能有交点
            return PASSBY;
    } else {// 两者均不为0且异号
        Point3d maybe_hP;
        double ratio = std::abs(dist_from) / (std::abs(dist_from) + std::abs(dist_to));
        maybe_hP = from + ratio * seg;
        if (testPointTriangle(maybe_hP, tri) == IN){
            hitPoint = maybe_hP;
            return HITTHROUGH;
        } else
            return PASSBY;
    }
}

TriangleTriangleType testTriangleTriangle(Triangle const& tri1, Triangle const& tri2, std::vector<Point3d> & hitPoints) {
    // 检查tri1的三条边是否与tri2相交
    SegmentTriangleType tri1_to_tri2[3];
    Point3d maybe_hp;
    for (int edge = 0; edge < 3; edge++) {
        tri1_to_tri2[edge] = testSegmentTriangle(tri1.m_pt[edge], tri1.m_pt[(edge+1)%3], tri2, maybe_hp);
        if (tri1_to_tri2[edge] == LYON) {
            hitPoints.clear();// 保证返回的hitPoints有且仅有2个点
            hitPoints.push_back(tri1.m_pt[edge]);
            hitPoints.push_back(tri1.m_pt[(edge+1)%3]);
            return INTERSECT;
        } else if (tri1_to_tri2[edge] == HITTHROUGH) {
            hitPoints.push_back(maybe_hp);
            if (hitPoints.size() == 2) 
                return INTERSECT;
        }// 线段与tri2为PASSBY不用处理
    }
    // 检查tri2的三条边是否与tri1相交
    SegmentTriangleType tri2_to_tri[3];
    for (int edge = 0; edge < 3; edge++) {
        tri2_to_tri[edge] = testSegmentTriangle(tri2.m_pt[edge], tri2.m_pt[(edge+1)%3], tri1, maybe_hp);
        if (tri2_to_tri[edge] == LYON) {
            hitPoints.clear();// 保证返回的hitPoints有且仅有2个点
            hitPoints.push_back(tri2.m_pt[edge]);
            hitPoints.push_back(tri2.m_pt[(edge+1)%3]);
            return INTERSECT;
        } else if (tri2_to_tri[edge] == HITTHROUGH) {
            hitPoints.push_back(maybe_hp);
            if (hitPoints.size() == 2) 
                return INTERSECT;
        }// 线段与tri2为PASSBY不用处理
    }
    assert(hitPoints.size() < 2);
    if (hitPoints.size() == 1) return TOUCH; 
    else return SEGREGATE;
}

