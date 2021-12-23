#ifndef __CAGEGENERATOR_H__
#define __CAGEGENERATOR_H__

#include <stdio.h>

#include <igl/readPLY.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/writePLY.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>

#include <Eigen/QR>
#include <Eigen/Sparse>
#include <igl/cotmatrix.h>
#include "igl/ray_mesh_intersect.h"
#include "igl/invert_diag.h"
#include <climits>
#include "Geometry.h"
//#include "igl/copyleft/cgal/intersect_other.h"

using namespace Eigen;
using namespace std;

class CageGenerator{
private:
    double sparseness_;
    MatrixXd V_mesh_,V_pca_basis_, V_cage_, V_cage_smooth_;
    MatrixXi F_mesh_, F_cage_;
    int num_vertices_mesh_, num_faces_mesh_;
    int num_vertices_cage_, num_faces_cage_;
    Vector3d barycenter_;// 所有vertex坐标的平均值，相当于mesh的中心
    MatrixXd pca_basis_matrix_;
    Vector3d bbBox_ptMin, bbBox_ptMax;// 特征空间下的点坐标
    int* n_vox;
    double* res_vox,* epsilon_tests;
    int num_feat_voxels;
    vector<vector<double> > voxel_center_pos_, voxel_grid_pos_;// voxel中心位置的坐标 voxel各个顶点的位置的坐标(均为特征空间) 紧致存储：高维为x,y,z方向，低维为该方向下标
    MatrixXd voxel_grid_pos_1Darray;
    MatrixXd voxel_face_normals_;
    vector < vector < vector<int> > > voxel_types_;
    vector < vector < vector<vec3d> > > voxel_grad_;

    MatrixXd outer_mesh_vertices_;
    MatrixXi outer_mesh_faces_;
    SparseMatrix<double> laplace_beltrami_matrix_;
    int num_smoothing_iterations_;
    double lambda_smooth_;
    
public:
    CageGenerator(const MatrixXd& vertices, const MatrixXi& faces, const int& smooth_it, const double& lambda,  const double& sparseness);
    void ComputePCA();
    void ComputePrincipalDirectionsBoundingBox(MatrixXd& bb_vertices, MatrixXi& bb_faces);
    void ComputeCubeMesh(MatrixXd& bb_vertices, MatrixXi& bb_faces, const Vector3d& min_point, const Vector3d& max_point);

    bool VoxelTriangleIntersection(const Triangle& triangle, const int i, const int j, const int k);
    void InitializeVoxelPositions();
    void InitializeVoxelNormals();
    int TupleToIndex(const int& i, const int& j, const int& k);

    int FindFeatureVoxels();
    int Modify_2Manifold();// 修正Course Bounding Cage满足2-manifold
    void ComputeSeparatingFaces(const Vector3i& curr_vox_tuple, const int& q,
                                Vector3i& face_1, Vector3i& face_2, const bool& reverse);
    void ExtractOuterSurface();
    void ComputeVoxelGrad();

    void SimplifyAdjacentFaces();
    void FillingAlgorithm(const Vector3i& vox_ind, vector<Vector3i>& outer_faces);
    void ComputeCage(MatrixXd& bb_vertices, MatrixXi& bb_faces);
    MatrixXd GetGridVertices();
    void SmoothCage();
    MatrixXd GetCage();
    void SetSmoothingParameters(const int& n_i, const float& lambda);
    MatrixXd GetSmoothCage();
    MatrixXi GetCageFaces();
};

class Wanted {
public:
    int i, j, k;
    int wanted_num;
    Wanted():i(0),j(0),k(0),wanted_num(0) {}
    Wanted(int ii, int jj, int kk, int wn):i(ii), j(jj), k(kk), wanted_num(wn) {}
    bool operator < (Wanted const & b) const {
        if (wanted_num == b.wanted_num) return i < b.i;
        return wanted_num > b.wanted_num;
    }
};

#endif