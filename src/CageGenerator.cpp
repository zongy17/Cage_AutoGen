#include "CageGenerator.h"

// 一个立方体有6个面
#define NUM_FACES_CUBE 6
// 一个三角形面片组成的立方体有12个面
#define NUM_FACES_CUBE_TRI_MESH 12
// 一个立方体有8个顶点
#define NUM_VERTICES_CUBE 8

static const double EPSILON_TESTS_RATIO = 0.001;

static bool implicit_smoothing = true;
static const double EPSILON_BOUNDING_BOX_ = 0.15;  // make the bounding box a bit larger

// For automatic cage generation : CHANGE THIS for a sparser or denser cage : the larger the parameter the denser the cage
CageGenerator::CageGenerator(const MatrixXd& vertices, const MatrixXi& faces, const int& smooth_it, \
                            const double& lambda, const double& sparseness) {
    sparseness_ = sparseness;
    V_mesh_ = vertices;
    F_mesh_ = faces;
    num_vertices_mesh_ = vertices.rows();
    num_faces_mesh_ = faces.rows();
    num_smoothing_iterations_ = smooth_it;
    // 计算这个mesh的所有vertex的坐标的算术平均
    barycenter_ = Vector3d::Zero(3);
    for (int i = 0; i < num_vertices_mesh_; i++){
        barycenter_ += vertices.row(i);
    }
    barycenter_ /= num_vertices_mesh_;
    InitializeVoxelNormals();
    lambda_smooth_ = lambda;
}

// 初始化一个体素的6个方向的法向量
void CageGenerator::InitializeVoxelNormals() {
    voxel_face_normals_ = MatrixXd::Zero(NUM_FACES_CUBE,3);
    voxel_face_normals_.row(0) = Vector3d::UnitX();
    voxel_face_normals_.row(1) = -Vector3d::UnitX();
    voxel_face_normals_.row(2) = Vector3d::UnitY();
    voxel_face_normals_.row(3) = -Vector3d::UnitY();
    voxel_face_normals_.row(4) = Vector3d::UnitZ();
    voxel_face_normals_.row(5) = -Vector3d::UnitZ();
}

// 根据bounding box在3个维度上的分辨率和顶点个数 生成各个voxel的中心位置grid_face_vertices_和它的顶点位置grid_face_vertices_
void CageGenerator::InitializeVoxelPositions() {
    assert(voxel_center_pos_.empty());
    assert(voxel_grid_pos_.empty());
    // 紧致存储
    for (int dir = 0; dir < 3; dir++) {
        voxel_center_pos_.push_back(vector<double>(n_vox[dir]));
        for (int idx = 0; idx < n_vox[dir]; idx++){
            double loc = bbBox_ptMin(dir) + ( ((double)idx + 0.5 )/ n_vox[dir] ) * (bbBox_ptMax(dir) - bbBox_ptMin(dir));
            voxel_center_pos_[dir][idx] = loc;
        }
    }
    for (int dir = 0; dir < 3; dir++) {
        voxel_grid_pos_.push_back(vector<double>(n_vox[dir]+1));
        for (int idx = 0; idx <= n_vox[dir]; idx++){
            double loc = bbBox_ptMin(dir) + ( ((double)idx)/ n_vox[dir] ) * (bbBox_ptMax(dir) - bbBox_ptMin(dir));
            voxel_grid_pos_[dir][idx] = loc;
        }
    }

    int num_vert_grid = (n_vox[0]+1) * (n_vox[1]+1) * (n_vox[2]+1);
    voxel_grid_pos_1Darray = MatrixXd::Zero(num_vert_grid, 3);
    for (int i = 0; i <= n_vox[0]; i++){// 第一个维度
        for (int j = 0; j <= n_vox[1]; j++){// 第二个维度
            for (int k = 0; k <= n_vox[2]; k++){// 第三个维度
                int ind = TupleToIndex(i, j, k);
                voxel_grid_pos_1Darray.row(ind) = Vector3d(voxel_grid_pos_[0][i], voxel_grid_pos_[1][j], voxel_grid_pos_[2][k]);
            }
        }
    }
    num_vertices_cage_ = num_vert_grid;// cage的顶点个数
}

// 根据特征空间中的min_point和max_point的坐标，以它俩为对角线构造一个cubeMesh
void CageGenerator::ComputeCubeMesh(MatrixXd &bb_vertices, MatrixXi &bb_faces, const Vector3d &min_point, const Vector3d &max_point) {
    int n_f = NUM_FACES_CUBE_TRI_MESH;
    bb_vertices = MatrixXd::Zero(NUM_VERTICES_CUBE, 3);// 最粗的bounding box顶点一共8个
    bb_faces = MatrixXi::Zero(n_f, 3);// 最粗的bounding box面一共12个（三角形面片）
    // Faces 初始化面的vertex id
    bb_faces.row(0) = Vector3i(0,1,2);
    bb_faces.row(1) = Vector3i(0,2,3);
    bb_faces.row(2) = Vector3i(1,5,6);
    bb_faces.row(3) = Vector3i(1,6,2);
    bb_faces.row(4) = Vector3i(5,4,7);
    bb_faces.row(5) = Vector3i(5,7,6);
    bb_faces.row(6) = Vector3i(4,0,3);
    bb_faces.row(7) = Vector3i(4,3,7);
    bb_faces.row(8) = Vector3i(4,5,1);
    bb_faces.row(9) = Vector3i(4,1,0);
    bb_faces.row(10) = Vector3i(3,2,6);
    bb_faces.row(11) = Vector3i(3,6,7);
    
    // Vertices
    //        7______________6
    //       /|             /|
    //     3/_|___________2/ |
    // z    | |           |  |
    // ^  y | |           |  |
    // |  ^ | |4__________|__|5
    // | /  | /           | /
    // |/   |/____________|/
    // |    0             1
    //  ------> x
    double xMin = min_point[0], yMin = min_point[1], zMin = min_point[2];
    double xMax = max_point[0], yMax = max_point[1], zMax = max_point[2];
    bb_vertices.row(0) = Vector3d(xMin, yMin, zMin);
    bb_vertices.row(1) = Vector3d(xMax, yMin, zMin);
    bb_vertices.row(2) = Vector3d(xMax, yMin, zMax);
    bb_vertices.row(3) = Vector3d(xMin, yMin, zMax);
    bb_vertices.row(4) = Vector3d(xMin, yMax, zMin);
    bb_vertices.row(5) = Vector3d(xMax, yMax, zMin);
    bb_vertices.row(6) = Vector3d(xMax, yMax, zMax);
    bb_vertices.row(7) = Vector3d(xMin, yMax, zMax);
}

// 计算mesh所有顶点构成的特征空间的基底，并将重心空间的坐标V_bar映射到特征空间去得到V_pca_basis_
void CageGenerator::ComputePCA() { 
    MatrixXd V_bar = V_mesh_;
    pca_basis_matrix_ = MatrixXd::Zero(3, 3);
    for (int i = 0; i < num_vertices_mesh_; i++){// V_bar是减去mesh中心得到的相对位移的坐标
        V_bar.row(i) -= barycenter_;
    }
    MatrixXd cov_mat = V_bar.transpose() * V_bar / num_vertices_mesh_;// V_bar顶点之间的协相关矩阵
    EigenSolver<MatrixXd> eigen_solver(cov_mat);// 用来求解协相关矩阵的特征向量
    // 特征方向提取
    for (int i = 0; i < 3; i++){
        Vector3d eigVec = eigen_solver.eigenvectors().col(i).real();
        pca_basis_matrix_.col(i) = eigVec / eigVec.norm();// 将单位化的特征向量作为PCA矩阵的方向向量
    }
    // pca_basis_matrix_的每一个列向量，都是特征空间的基向量，在当前重心空间（即以e0=(1,0,0),e1=(0,1,0),e2=(0,0,1)为基底）下的表示（当前基底的线性组合）
    // 即 [p0,p1,p2] = [e0,e1,e2] * pca_basis_matrix_
    // 在重心空间中的坐标 V_bar 和在特征空间中的坐标 V_pca_basis_ 应满足：[e0,e1,e2] * V_bar == [p0,p1,p2] * V_pca_basis_
    // 因此 V_pca_basis_ = [p0,p1,p2]^(-1) * [e0,e1,e2] * V_bar = pca_basis_matrix_^(-1) * V_bar
    // 而注意pca_basis_matrix_是一个正交阵，它的逆就是自身的转置
    // 即将重心空间中的顶点坐标V_bar还原出特征空间下的顶点坐标V_pca_basis
    V_pca_basis_ = (pca_basis_matrix_.transpose() * V_bar.transpose()).transpose();// 用PCA矩阵作用于V_bar
}

// 计算最大（最粗）的bounding box，根据PCA结果来确定box的三个维度的方向，返回在一般世界空间中的bb_vertices和bb_faces
void CageGenerator::ComputePrincipalDirectionsBoundingBox(MatrixXd &bb_vertices, MatrixXi &bb_faces) {
    ComputePCA();
    // Change basis
    res_vox = new double[3];
    n_vox = new int[3];
    // 找到PCA变换后最极端的顶点（离中心最远）
    double minCoords[3] = {__DBL_MAX__, __DBL_MAX__, __DBL_MAX__};
    double maxCoords[3] = {__DBL_MIN__, __DBL_MIN__, __DBL_MIN__};// 这里是不是搞错了？是不是应该用负值？-__DBL_MAX__
    // double maxCoords[3] = {-__DBL_MAX__, -__DBL_MAX__, -__DBL_MAX__};
    Vector3d point;
    for (int i = 0; i < num_vertices_mesh_; i++){
        point = V_pca_basis_.row(i);// V_pca_basis_以pca（特征方向）作为基底的坐标
        for (int j = 0; j < 3; j++){
            if (point(j) > maxCoords[j]){
                maxCoords[j] = point(j);
            }
            if (point(j) < minCoords[j]){
                minCoords[j] = point(j);
            }
        }
    }
    assert(minCoords[0]<0 && minCoords[1]<0 && minCoords[2]<0);
    assert(maxCoords[0]>0 && maxCoords[1]>0 && maxCoords[2]>0);

    bbBox_ptMin = Vector3d::Zero(3);
    bbBox_ptMax = Vector3d::Zero(3);
    epsilon_tests = new double[3];
    // Bounding Box在3维上的宽度（特征空间）
    Vector3d deltaBB(maxCoords[0] - minCoords[0], maxCoords[1] - minCoords[1], maxCoords[2] - minCoords[2]);
    // 找到bounding box宽度最宽的那一维
    int j_max;
    double dim_max = __DBL_MIN__;
    for (int j = 0; j < 3; j++){
        bbBox_ptMin(j) = minCoords[j] - EPSILON_BOUNDING_BOX_ * deltaBB(j);// 使得bounding box比mesh稍微大一点
        bbBox_ptMax(j) = maxCoords[j] + EPSILON_BOUNDING_BOX_ * deltaBB(j);
        double dim_size = bbBox_ptMax(j) - bbBox_ptMin(j);
        if (dim_size > dim_max){
            dim_max = dim_size;
            j_max = j;
        }
    }
    // 最宽的那一维的体素数目 number of voxels，根据稀疏程度计算
    n_vox[j_max] = (int)(sqrt(sparseness_ * double(num_vertices_mesh_) / 6.)); // Chuhua Xian. https://www.cse.wustl.edu/~taoju/research/meanvalue.pdf
    std::cout << "num_vox = " << n_vox[j_max] << " in " << j_max << "-th dimension" << std::endl;
    // 最宽的那一维的体素的分辨率 resolution of voxelization
    res_vox[j_max] = (bbBox_ptMax(j_max) - bbBox_ptMin(j_max)) / n_vox[j_max];
    for (int j = 0; j < 3; j++){
        if (j != j_max){// 余下两个非最宽的维度的体素的数目和分辨率
            int n = (int)((bbBox_ptMax(j) - bbBox_ptMin(j)) / res_vox[j_max]) + 2;
            n_vox[j] = n;
            std::cout << "num_vox = " << n_vox[j] << " in " << j << "-th dimension" << std::endl;
            res_vox[j] = (bbBox_ptMax(j) - bbBox_ptMin(j)) / n_vox[j];
            epsilon_tests[j] = EPSILON_TESTS_RATIO * res_vox[j];
        }
    }
    // 根据bbBox_ptMin, bbBox_ptMax （bounding box的最远点）计算特征空间下的bounding box的立方体mesh，顶点bb_vertices，面bb_faces
    ComputeCubeMesh(bb_vertices, bb_faces, bbBox_ptMin, bbBox_ptMax);
    // Change basis （PCA变换，还原回重心空间下的bounding box）
    bb_vertices = (pca_basis_matrix_ * bb_vertices.transpose()).transpose();
    // 加上mesh的中心，还原回一般世界空间的bounding box
    for (int i = 0; i < bb_vertices.rows(); i++){
        bb_vertices.row(i) = barycenter_ +  Vector3d(bb_vertices.row(i));
    }
}

// 检测原mesh的一个三角形面片mesh_tri(特征空间)，是否与(i,j,k)处的voxel有交集
bool CageGenerator::VoxelTriangleIntersection(const Triangle& mesh_tri, const int i, const int j, const int k) {
    // Vertices
    //       (011)
    //        7______________6(111)
    // (001) /|      (101)  /|
    //     3/_|___________2/ |
    // z    | |           |  |
    // ^  y | | (010)     |  |
    // |  ^ | |4__________|__|5(110)
    // | /  | /           | /
    // |/   |/____________|/
    //      0(000)         1(100)
    //  ------> x
    static const int ngb_idx_offset[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1},
        {0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {0, 1, 1}
    };
    static const int face_ngb_indices[12][2] = {
        {0, 1}, {1, 5}, {5, 4}, {4, 0},// k相同的
        {3, 2}, {2, 6}, {6, 7}, {7, 3},// k+1相同的
        {0, 3}, {1, 2}, {5, 6}, {4, 7} // i和j相同的
    };

    Point3d center(voxel_center_pos_[0][i], voxel_center_pos_[1][j], voxel_center_pos_[2][k]);
    for (unsigned f = 0; f < 12; f++) {
        // printf("f %d\n", f);
        int f_ngb1_idx = face_ngb_indices[f][0];
        int f_ngb2_idx = face_ngb_indices[f][1];
        int p1_i = i + ngb_idx_offset[f_ngb1_idx][0], p1_j = j + ngb_idx_offset[f_ngb1_idx][1], p1_k = k + ngb_idx_offset[f_ngb1_idx][2];
        int p2_i = i + ngb_idx_offset[f_ngb2_idx][0], p2_j = j + ngb_idx_offset[f_ngb2_idx][1], p2_k = k + ngb_idx_offset[f_ngb2_idx][2];

        Point3d p_grid_1(voxel_grid_pos_[0][p1_i], voxel_grid_pos_[1][p1_j], voxel_grid_pos_[2][p1_k]);
        Point3d p_grid_2(voxel_grid_pos_[0][p2_i], voxel_grid_pos_[1][p2_j], voxel_grid_pos_[2][p2_k]);
        Triangle voxel_tri(center, p_grid_1, p_grid_2);
        std::vector<Point3d> interSect_pts; interSect_pts.clear();
        // double max_norm_dist = __DBL_MAX__;
        // if (TriangleIntersectionTest(voxel_tri, mesh_tri, interSect_pts) == INTERSECTION) {
        //     // assert(interSect_pts.size() == 2);
        //     if (interSect_pts.size() != 2) {
        //         printf("i %d j %d k %d center %lf %lf %lf\n",    i,    j,    k, center(0), center(1), center(2));
        //         printf("i %d j %d k %d pgrid1 %lf %lf %lf\n", p1_i, p1_j, p1_k, p_grid_1(0), p_grid_1(1), p_grid_1(2));
        //         printf("i %d j %d k %d pgrid2 %lf %lf %lf\n", p2_i, p2_j, p2_k, p_grid_2(0), p_grid_2(1), p_grid_2(2));
        //         printf(" mesh tri 0: %lf %lf %lf\n", mesh_tri.m_pt[0](0), mesh_tri.m_pt[0](1), mesh_tri.m_pt[0](2));
        //         printf(" mesh tri 1: %lf %lf %lf\n", mesh_tri.m_pt[1](0), mesh_tri.m_pt[1](1), mesh_tri.m_pt[1](2));
        //         printf(" mesh tri 2: %lf %lf %lf\n", mesh_tri.m_pt[2](0), mesh_tri.m_pt[2](1), mesh_tri.m_pt[2](2));
        //         printf(" size: %d\n", interSect_pts.size());
        //         for (int it = 0; it < interSect_pts.size(); it++)
        //             printf(" %lf %lf %lf\n", interSect_pts[it](0), interSect_pts[it](1), interSect_pts[it](2));
        //         std::vector<Point3d> test; test.clear();
        //         TriangleIntersectionTest(voxel_tri, mesh_tri, test);
        //         exit(1);
        //     }
        //     max_norm_dist = std::min(max_norm_dist, (center-mesh_tri.m_pt[0]).norm());// O - P0
        //     max_norm_dist = std::min(max_norm_dist, (center-mesh_tri.m_pt[1]).norm());// O - P1
        //     max_norm_dist = std::min(max_norm_dist, (center-mesh_tri.m_pt[2]).norm());// O - P2
        //     max_norm_dist = std::min(max_norm_dist, (center-interSect_pts[0]).norm());// O - t1
        //     max_norm_dist = std::min(max_norm_dist, (center-interSect_pts[1]).norm());// O - t2
        //     double w_v;
        //     if      (p1_i == p2_i) w_v = res_vox[0];
        //     else if (p1_j == p2_j) w_v = res_vox[1];
        //     else if (p1_k == p2_k) w_v = res_vox[2];
        //     else {
        //         printf(" p1: i %d j %d k %d\n p2: i %d j %d k %d\n", p1_i, p1_j, p1_k, p2_i, p2_j, p2_k);
        //         exit(1);
        //     }
        //     if (max_norm_dist < 0.5*w_v) return true;
        // }
        if (testTriangleTriangle(mesh_tri, voxel_tri, interSect_pts) != SEGREGATE) {
            // printf("i %d j %d k %d center %lf %lf %lf\n",    i,    j,    k, center(0), center(1), center(2));
            // printf("i %d j %d k %d pgrid1 %lf %lf %lf\n", p1_i, p1_j, p1_k, p_grid_1(0), p_grid_1(1), p_grid_1(2));
            // printf("i %d j %d k %d pgrid2 %lf %lf %lf\n", p2_i, p2_j, p2_k, p_grid_2(0), p_grid_2(1), p_grid_2(2));
            // printf(" mesh tri 0: %lf %lf %lf\n", mesh_tri.m_pt[0](0), mesh_tri.m_pt[0](1), mesh_tri.m_pt[0](2));
            // printf(" mesh tri 1: %lf %lf %lf\n", mesh_tri.m_pt[1](0), mesh_tri.m_pt[1](1), mesh_tri.m_pt[1](2));
            // printf(" mesh tri 2: %lf %lf %lf\n", mesh_tri.m_pt[2](0), mesh_tri.m_pt[2](1), mesh_tri.m_pt[2](2));
            return true;
        }
    }
    return false;
}

// the voxels intersecting the mesh M are taken as the feature voxels: voxel_types_[i][j][k] = 0
// the other voxels are categorized as inner voxels: voxel_types_[i][j][k] = -1 at the beginning
// 每个feature voxel都生成一个cubeMesh，并将其顶点和面存储于feat_voxels_vertices和feat_voxels_faces，注意面的朝向为反的
int CageGenerator::FindFeatureVoxels() { 
    assert(voxel_types_.empty());
    voxel_types_ = vector<vector<vector<int>>>(n_vox[0], vector<vector<int>>(n_vox[1], vector<int>(n_vox[2], -1)));// 初始化为-1: inner

    int cnt_feat_voxel = 0;
    Vector3d deltaVox(res_vox[0]/2., res_vox[1]/2., res_vox[2]/2.);
    // 换一种思路，对mesh的所有面face进行遍历，定位face附近的voxel，然后逐个“注册”为feature voxel
    for (int f = 0; f < num_faces_mesh_; f++) {
        Triangle mesh_tri(V_pca_basis_.row(F_mesh_(f,0)), V_pca_basis_.row(F_mesh_(f,1)), V_pca_basis_.row(F_mesh_(f,2)));
        // 确定这个三角形面片在三个维度上最小最大能去到什么位置
        int search_beg[3], search_end[3];
        double minCoords[3] = {mesh_tri.m_pt[0].x(), mesh_tri.m_pt[0].y(), mesh_tri.m_pt[0].z()};
        double maxCoords[3] = {mesh_tri.m_pt[0].x(), mesh_tri.m_pt[0].y(), mesh_tri.m_pt[0].z()};
        for (int dir = 0; dir < 3; dir++){
            minCoords[dir] = std::min(minCoords[dir], std::min(mesh_tri.m_pt[1](dir), mesh_tri.m_pt[2](dir)));
            maxCoords[dir] = std::max(minCoords[dir], std::max(mesh_tri.m_pt[1](dir), mesh_tri.m_pt[2](dir)));
            // upper_bound返回第一个大于val的迭代器，而lower_bound是第一个大于等于
            search_beg[dir] = (upper_bound(voxel_grid_pos_[dir].begin(), voxel_grid_pos_[dir].end(), minCoords[dir]) - voxel_grid_pos_[dir].begin()) - 1;
            assert(search_beg[dir] >= 0 && search_beg[dir] < n_vox[dir]);
            search_end[dir] = (lower_bound(voxel_grid_pos_[dir].begin(), voxel_grid_pos_[dir].end(), maxCoords[dir]) - voxel_grid_pos_[dir].begin()) - 1;
            assert(search_end[dir] >= 0 && search_end[dir] < n_vox[dir]);
        }
        // 根据是否相交判断是否为feature voxel
        for (int i = search_beg[0]; i <= search_end[0]; i++) {// 注意是 i<=
            for (int j = search_beg[1]; j <= search_end[1]; j++) {
                for (int k = search_beg[2]; k <= search_end[2]; k++) {
                    int & this_voxel_type = voxel_types_[i][j][k];
                    if (this_voxel_type == 0) continue;// 已经“注册”过了,不用再检测了
                    
                    if (VoxelTriangleIntersection(mesh_tri, i, j, k)) {// 检测这个三角形面片是否与voxel有交集
                        this_voxel_type = 0;
                        cnt_feat_voxel++;
                    }
                }// k
            }// j
        }// i
    }// faces

    // 修补CBC满足2-manifold
    cnt_feat_voxel += Modify_2Manifold();

    std::cout << "num_feature_voxels = " << cnt_feat_voxel << std::endl;
    return cnt_feat_voxel;
}

// 获取以curr_vox_tuple为下标的voxel在它的q方向（6个方向之一）的邻居voxel之间的分割面
// 如果reverse=false，则返回的两个面face_1, face_2的朝向为邻居voxel指向curr_voxel
void CageGenerator::ComputeSeparatingFaces(const Vector3i& curr_vox_tuple, const int& q,
                                           Vector3i& face_1, Vector3i& face_2, const bool& reverse) {
    int i = curr_vox_tuple(0), j = curr_vox_tuple(1), k = curr_vox_tuple(2);
    int i0, i1, i2, i3;
    int j0, j1, j2, j3;
    int k0, k1, k2, k3;
    // Separating on X-axis
    if(q == 0){//         2
        //               /|
        //             3/ |
        //              | | 
        // curr_voxel   | |  neighbor_voxel
        //   i(full)    | |   i+1(full)
        // ^z  ^y       | /0
        // |  /        1|/
        // | /        i+1(half)
        // -------------------------> x
        i0 = i1 = i2 = i3 = i+1;
        k0 = k1 = k;
        k2 = k3 = k + 1;
        j1 = j3 = j;
        j0 = j2 = j + 1;
    }
    else if (q==1){//         3
        //                   /|
        //                 2/ |
        //                  | | 
        // neighbor_voxel   | |  curr_voxel
        //     i-1(full)    | |   i (full)
        // ^z  ^y           | /1
        // |  /            0|/
        // | /            i(half)
        // -------------------------> x

        i0 = i1 = i2 = i3 = i;
        k0 = k1 = k;
        k2 = k3 = k + 1;
        j1 = j3 = j + 1;
        j0 = j2 = j;
    }
    
    // Separating on Y-axis
    if(q == 2){//         3
        //               /|
        //             1/ |
        //              | | 
        // curr_voxel   | |  neighbor_voxel
        //   j(full)    | |   j+1(full)
        // ^x  ^z       | /2
        // |  /        0|/
        // | /        j+1(half)
        // -------------------------> y
        j0 = j1 = j2 = j3 = j+1;
        i0 = i2 = i;
        i1 = i3 = i+1;
        k0 = k1 = k;
        k2 = k3 = k+1;
    }
    else if (q==3){//         2
        //                   /|
        //                 0/ |
        //                  | | 
        // neighbor_voxel   | |  curr_voxel
        //     j-1(full)    | |   j (full)
        // ^x  ^z           | /3
        // |  /            1|/
        // | /            j(half)
        // -------------------------> y
        j0 = j1 = j2 = j3 = j;
        i0 = i2 = i+1;
        i1 = i3 = i;
        k0 = k1 = k;
        k2 = k3 = k+1;
    }
    
    // Separating on Z-axis
    if(q == 4){//         1
        //               /|
        //             0/ |
        //              | | 
        // curr_voxel   | |  neighbor_voxel
        //   k(full)    | |   k+1(full)
        // ^y  ^x       | /3
        // |  /        2|/
        // | /        k+1(half)
        // -------------------------> z
        k0 = k1 = k2 = k3 = k+1;
        i0 = i2 = i;
        i1 = i3 = i+1;
        j2 = j3 = j;
        j0 = j1 = j+1;
    }
    else if (q==5){//         0
        //                   /|
        //                 1/ |
        //                  | | 
        // neighbor_voxel   | |  curr_voxel
        //     k-1(full)    | |   k (full)
        // ^y  ^x           | /2
        // |  /            3|/
        // | /            k(half)
        // -------------------------> z
        k0 = k1 = k2 = k3 = k;
        i0 = i2 = i+1;
        i1 = i3 = i;
        j2 = j3 = j;
        j0 = j1 = j+1;
    }
    
    face_1 = Vector3i::Zero(3), face_2 = Vector3i::Zero(3);
    int ind1 = 1, ind2 = 2;
    if (reverse){
        ind1 = 2;
        ind2 = 1;
    }
    face_1(0) = TupleToIndex(i0,j0,k0);
    face_1(ind1) = TupleToIndex(i1,j1,k1);
    face_1(ind2) = TupleToIndex(i2,j2,k2);
    face_2(0) = TupleToIndex(i1,j1,k1);
    face_2(ind1) = TupleToIndex(i3,j3,k3);
    face_2(ind2) = TupleToIndex(i2,j2,k2);
}

// 提取外向面
void CageGenerator::ExtractOuterSurface() {
    int i,j,q;
    Vector3i seed_index_tuple;
    vector<Vector3i> outer_faces;
    outer_faces.clear();
    // Expand from seeds 从bounding box的六个顶点上的voxel开始深度优先搜索
    for (int k = 0; k < 2; k++){
        i = k * (n_vox[0] - 1);
        for (int l = 0; l < 2; l++){
            j = l * (n_vox[1] - 1);
            for (int r = 0; r < 2; r++){
                q = r * (n_vox[2] - 1);
                seed_index_tuple = Vector3i(i,j,q);
                FillingAlgorithm(seed_index_tuple, outer_faces);
            }
        }
    }
    // std::cout << "Done filling from boundary." << std::endl;

    // 做完从6个顶点开始的Filling之后，应该有：
    // mesh内部的inner voxel仍然是-1
    // feature voxel仍然是0
    // mesh外部的outer voxel标记为1
    // 此时outer_faces里的分割面的朝向均为从自己指向邻居voxel(feature voxel)

    // 接下去处理位于Bounding Box边界上的outer_faces（之前的FillingAlgorithm无法覆盖这种情况），原则也是分割面的朝向指向feature voxel
    int voxel_info;
    Vector3i curr_vox_tuple, neighbor_vox_tuple, face_1, face_2;
    bool is_neighbor_inside_bb;
    for (int i = 0; i < n_vox[0]; i++)
        for (int j = 0; j < n_vox[1]; j++)
            for (int k = 0; k < n_vox[2]; k++){
                voxel_info = voxel_types_[i][j][k];
                curr_vox_tuple = Vector3i(i,j,k);
                if (voxel_info == 0){// feature voxel
                    for (int q = 0; q < NUM_FACES_CUBE; q++){// 遍历这个feature voxel在6个方向上的邻居
                        neighbor_vox_tuple = curr_vox_tuple + Vector3i(voxel_face_normals_.row(q).cast<int>());// 方向q的邻居voxel的三元id
                        is_neighbor_inside_bb = true;
                        for (int l = 0; l < 3; l++){// 检查是否越界
                            if (neighbor_vox_tuple(l) < 0 || neighbor_vox_tuple(l) > n_vox[l] - 1){
                                is_neighbor_inside_bb = false;// 如果越界说明这个voxel自身处于Bounding Box的边界面上
                                break;
                            }
                        }
                        if (!is_neighbor_inside_bb){
                            // 则也要创建分割面，并且此时分割面的朝向指向自己
                            ComputeSeparatingFaces(curr_vox_tuple,q, face_1, face_2, false);
                            outer_faces.push_back(face_1);
                            outer_faces.push_back(face_2);
                        }
                    }
                }
            }

    // 根据outer_faces的信息整合成cage
    num_faces_cage_ = outer_faces.size();// cage的面就是outer_faces
    std::cout << "Number of outward faces after cage extraction : " << num_faces_cage_ << std::endl;
    F_cage_ = MatrixXi::Zero(num_faces_cage_,3);
    
    // 根据外向面构造“规规整整”的三角形网格：Construct mesh based on outer faces
    vector<Vector3d> vector_vertices;
    vector_vertices.clear();
    VectorXi added_vertices = - VectorXi::Ones(voxel_grid_pos_1Darray.rows());// grid_face_vertices顶点先全部置为-1，谁被加入了，就改为加入的序号
    int count = 0;// 加入的序号/当前总数
    for (int f = 0; f < num_faces_cage_; f++){// 遍历每个outer_face
        Vector3i face = outer_faces[f];// 组成face的顶点idx是用voxelization得到的grid的顶点来表示的
        Vector3i new_face = Vector3i::Zero(3);
        for (int i = 0; i < 3; i++){
            int ind = face(i);// outer_face的一个顶点a的index
            int curr_vert_ind = added_vertices(ind);// 如果这个顶点a还没有被加入
            if (curr_vert_ind == -1){
                added_vertices(ind) = count;
                new_face(i) = count;
                vector_vertices.push_back(voxel_grid_pos_1Darray.row(ind));// 记录这个顶点a的位置坐标
                count++;
            }
            else{// 如果这个顶点a已经被加入过了
                new_face(i) = curr_vert_ind;
            }
        }
        F_cage_.row(f) = new_face;
    }
    // 以上的construction相当于只是将outer_Faces换了一种表达方式：
    // 用voxelization生成的三维规则格点序号（三维i,j,k <=>一维线性）表示，转换成仅用outer_faces所含的顶点来表示
    // 相当于其余outer voxel和inner voxel的信息都可以丢弃掉了
    num_vertices_cage_ = count;
    V_cage_ = MatrixXd::Zero(num_vertices_cage_,3);
    for (int i = 0; i < num_vertices_cage_; i++){
        V_cage_.row(i) = vector_vertices.at(i);
    }
}

void CageGenerator::FillingAlgorithm(const Vector3i &curr_vox_tuple, vector<Vector3i> &outer_faces) {
    Vector3i neighbor_vox_tuple;
    int i = curr_vox_tuple(0), j = curr_vox_tuple(1), k = curr_vox_tuple(2);
    int i_neighbor, j_neighbor, k_neighbor;
    int neighbor_partition_info;
    bool is_neighbor_inside_bb;
    Vector3i face_1, face_2;
    // Mark as visited
    voxel_types_[i][j][k] = 1;// outer voxels
    for (int q = 0; q < NUM_FACES_CUBE; q++){// 遍历这个voxel在6个方向上的邻居
        neighbor_vox_tuple = curr_vox_tuple + Vector3i(voxel_face_normals_.row(q).cast<int>());
        // Check that neighboring voxel is inside the BB
        is_neighbor_inside_bb = true;
        for (int l = 0; l < 3; l++){// 检查是否越界
            if (neighbor_vox_tuple(l) < 0 || neighbor_vox_tuple(l) > n_vox[l] - 1){
                is_neighbor_inside_bb = false;
                break;
            }
        }// break for loop
        if (is_neighbor_inside_bb){
            i_neighbor = neighbor_vox_tuple(0);
            j_neighbor = neighbor_vox_tuple(1);
            k_neighbor = neighbor_vox_tuple(2);
            neighbor_partition_info = voxel_types_[i_neighbor][j_neighbor][k_neighbor];
            if (neighbor_partition_info == 1){
                // already visited
                continue;// continue for (int q; ...) loop
            } else if(neighbor_partition_info == 0){// 这个voxel的邻居voxel是一个feature voxel！
                // Feature voxel : Add the separating face
                // 计算自己（这个voxel）和邻居（feature voxel）的分割面（亦即两个voxel的公共面）
                ComputeSeparatingFaces(curr_vox_tuple,q, face_1, face_2, true);
                outer_faces.push_back(face_1);
                outer_faces.push_back(face_2);
            } else{// neighbor_partition_info == -1
                // Visit neighbor: Depth-First Search!
                FillingAlgorithm(neighbor_vox_tuple, outer_faces);
            }
        }
    }
}

int CageGenerator::Modify_2Manifold() {
    int cnt_mod_feat_voxel = 0;
    // 修正non-manifold edges：Voxel-attaching operation
    for (int i = 0; i < n_vox[0]; i++)
        for (int j = 0; j < n_vox[1]; j++)
            for (int k = 0; k < n_vox[2]; k++){
                if (voxel_types_[i][j][k] != 0) continue;// 非feature voxel，不用管

                bool non_2_manifold[3][3][3] = {{false, false, false}, {false, false, false}, {false, false, false}};
                int wanted[3][3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
                // x 不动
                for (int jo = -1; jo <= 1; jo += 2)
                    for (int ko = -1; ko <= 1; ko += 2){
                        int ngb_j = j + jo;
                        int ngb_k = k + ko;
                        if (ngb_j>=0 && ngb_j<n_vox[1] && ngb_k>=0 && ngb_k<n_vox[2]) {// 二维斜角邻居在整个voxel网格内
                            if (voxel_types_[i][ngb_j][ngb_k] != 0) continue;// 二维斜角邻居不是feature voxel，不用管
                            // assert(voxel_types_[i][ngb_j][ngb_k]==0);
                            if (voxel_types_[i][j][ngb_k] != 0 && voxel_types_[i][ngb_j][k] != 0) {
                                non_2_manifold[0+1][jo+1][ko+1] = true;
                                wanted[0+1][0+1][ko+1]++;
                                wanted[0+1][jo+1][0+1]++;
                            }
                        }
                    }
                // y 不动
                for (int io = -1; io <= 1; io += 2)
                    for (int ko = -1; ko <= 1; ko += 2){
                        int ngb_i = i + io;
                        int ngb_k = k + ko;
                        if (ngb_i>=0 && ngb_i<n_vox[0] && ngb_k>=0 && ngb_k<n_vox[2]) {// 二维斜角邻居在整个voxel网格内
                            if (voxel_types_[ngb_i][j][ngb_k] != 0) continue;// 二维斜角邻居不是feature voxel，不用管
                            // assert(voxel_types_[ngb_i][j][ngb_k]==0);
                            if (voxel_types_[i][j][ngb_k] != 0 && voxel_types_[ngb_i][j][k] != 0) {
                                non_2_manifold[io+1][0+1][ko+1] = true;
                                wanted[0+1][0+1][ko+1]++;
                                wanted[io+1][0+1][0+1]++;
                            }
                        }
                    }
                // z 不动
                for (int io = -1; io <= 1; io += 2)
                    for (int jo = -1; jo <= 1; jo += 2){
                        int ngb_i = i + io;
                        int ngb_j = j + jo;
                        if (ngb_i>=0 && ngb_i<n_vox[0] && ngb_j>=0 && ngb_j<n_vox[1]) {// 二维斜角邻居在整个voxel网格内
                            if (voxel_types_[ngb_i][ngb_j][k] != 0) continue;// 二维斜角邻居不是feature voxel，不用管
                            // assert(voxel_types_[ngb_i][ngb_j][k]==0);
                            if (voxel_types_[ngb_i][j][k] != 0 && voxel_types_[i][ngb_j][k] != 0) {
                                non_2_manifold[io+1][jo+1][0+1] = true;
                                wanted[0+1][jo+1][0+1]++;
                                wanted[io+1][0+1][0+1]++;
                            }
                        }
                    }
                std::vector<Wanted> wanted_arr = {\
                    Wanted(-1, 0, 0, wanted[0][1][1]), Wanted(1, 0, 0, wanted[2][1][1]),\
                    Wanted( 0,-1, 0, wanted[1][0][1]), Wanted(0, 1, 0, wanted[1][2][1]),\
                    Wanted( 0, 0,-1, wanted[1][1][0]), Wanted(0, 0, 1, wanted[1][1][2]) \
                };
                assert(wanted_arr.size()==6);
                std::sort(wanted_arr.begin(), wanted_arr.end());// 升序排列
                int ptr = 0;
                while (ptr < 6 && wanted_arr[ptr].wanted_num > 0) {
                    int io = wanted_arr[ptr].i;
                    int jo = wanted_arr[ptr].j;
                    int ko = wanted_arr[ptr].k;
                    bool need_to_add = false;
                    if (io != 0) {
                        assert(jo==0 && ko==0);
                        if (non_2_manifold[io+1][1][0] || non_2_manifold[io+1][1][2] || non_2_manifold[io+1][0][1] || non_2_manifold[io+1][2][1]) {
                            non_2_manifold[io+1][1][0] =  non_2_manifold[io+1][1][2] =  non_2_manifold[io+1][0][1] =  non_2_manifold[io+1][2][1] = false;
                            need_to_add = true;
                        }
                    } else if (jo != 0) {
                        assert(io==0 && ko==0);
                        if (non_2_manifold[1][jo+1][0] || non_2_manifold[1][jo+1][2] || non_2_manifold[0][jo+1][1] || non_2_manifold[2][jo+1][1]) {
                            non_2_manifold[1][jo+1][0] =  non_2_manifold[1][jo+1][2] =  non_2_manifold[0][jo+1][1] =  non_2_manifold[2][jo+1][1] = false;
                            need_to_add = true;
                        }
                    } else if (ko != 0) {
                        assert(io==0 && jo==0);
                        if (non_2_manifold[1][0][ko+1] || non_2_manifold[1][2][ko+1] || non_2_manifold[0][1][ko+1] || non_2_manifold[2][1][ko+1]) {
                            non_2_manifold[1][0][ko+1] =  non_2_manifold[1][2][ko+1] =  non_2_manifold[0][1][ko+1] =  non_2_manifold[2][1][ko+1] = false;
                            need_to_add = true;
                        }
                    } else {
                        printf("Error: io %d jo %d ko %d wanted_num %d\n", io, jo, ko, wanted_arr[ptr].wanted_num);
                        exit(1);
                    }
                    // 添加这个位置pos作为feature voxel
                    if (need_to_add) {
                        voxel_types_[i+io][j+jo][k+ko] = 0;
                        cnt_mod_feat_voxel++;
                    }
                    ptr++;
                }
                // 做完之后应该保证没有non_manifold了
                assert((!non_2_manifold[1][0][0]) && (!non_2_manifold[1][0][2]) && (!non_2_manifold[1][2][0]) && (!non_2_manifold[1][2][2]) \
                     &&(!non_2_manifold[0][1][0]) && (!non_2_manifold[0][1][2]) && (!non_2_manifold[2][1][0]) && (!non_2_manifold[2][1][2]) \
                     &&(!non_2_manifold[0][0][1]) && (!non_2_manifold[0][2][1]) && (!non_2_manifold[2][0][1]) && (!non_2_manifold[2][2][1]) );
            }
    return cnt_mod_feat_voxel;
}

void CageGenerator::ComputeVoxelGrad() {
    assert(voxel_grad_.empty());
    voxel_grad_ = vector<vector<vector<vec3d>>>(n_vox[0], vector<vector<vec3d>>(n_vox[1], vector<vec3d>(n_vox[2], vec3d(0.0, 0.0, 0.0))));// 初始化为0.0

    // 注意voxel_grad是定义在voxel体内的（full）
    for (int i = 0; i < n_vox[0]-1; i++) 
        for (int j = 0; j < n_vox[1]-1; j++)
            for (int k = 0; k < n_vox[2]-1; k++) {
                voxel_grad_[i][j][k] = vec3d(voxel_types_[i+1][j][k] - voxel_types_[i][j][k], \
                                             voxel_types_[i][j+1][k] - voxel_types_[i][k][k], \
                                             voxel_types_[i][j][k+1] - voxel_types_[i][j][k]);
            }
}

void CageGenerator::ComputeCage(MatrixXd& bb_vertices, MatrixXi& bb_faces) {
    ComputePrincipalDirectionsBoundingBox(bb_vertices, bb_faces);
    InitializeVoxelPositions();
    num_feat_voxels = FindFeatureVoxels();
    ExtractOuterSurface();
    igl::writeOBJ("cage_vox.obj", GetCage(), GetCageFaces());
    // ComputeVoxelGrad();
    //SimplifyAdjacentFaces();
    SmoothCage();
}

// 返回位于三维数组中下标为(i,j,k)的数据转换成一维后的线性下标：i + (nx+1)*(j + (ny+1)*k)
int CageGenerator::TupleToIndex(const int& i, const int& j, const int& k) {
    int res = k * (n_vox[0]+1) * (n_vox[1]+1) + j * (n_vox[0]+1) + i;
    return res;
}

// 返回cage的所有网格顶点（数目为(nvox[0]+1)*(nvox[1]+1)*(nvox[2]+1)）在世界空间中的坐标
MatrixXd CageGenerator::GetGridVertices() {
    // MatrixXd res = (pca_basis_matrix_ * grid_face_vertices_.transpose()).transpose();
    MatrixXd res = (pca_basis_matrix_ * voxel_grid_pos_1Darray.transpose()).transpose();
    for (int i = 0; i < res.rows(); i++){
        res.row(i) = barycenter_ + Vector3d(res.row(i));
    }
    return res;
}

// 返回cage的面对应的顶点索引
MatrixXi CageGenerator::GetCageFaces() { 
    return F_cage_;
}

void CageGenerator::SetSmoothingParameters(const int& n_i, const float& lambda) {
    lambda_smooth_ = lambda;
    num_smoothing_iterations_ = n_i;
}

void CageGenerator::SmoothCage() { 
    // Alternative discrete mean curvature
    // MatrixXd HN;
    // SparseMatrix<double> L,M,Minv, HN;
    // igl::cotmatrix(V_cage_, F_cage_, L);
    // igl::massmatrix(V_cage_, F_cage_, igl::MASSMATRIX_TYPE_VORONOI, M);
    // igl::invert_diag(M, Minv);
    // // Laplace-Beltrami of position
    // HN = -Minv*(L * V_cage_);
    // // Extract magnitude as mean curvature
    // VectorXd H = HN.rowwise().norm();
    
    igl::cotmatrix(V_cage_, F_cage_, laplace_beltrami_matrix_);
    // SparseMatrix<double> M, invM;
    // igl::massmatrix(V_cage_, F_cage_, igl::MASSMATRIX_TYPE_VORONOI, M);
    // igl::invert_diag(M, invM);
    // laplace_beltrami_matrix_ = invM * laplace_beltrami_matrix_;

    MatrixXd V_cage_prev = V_cage_;
    V_cage_smooth_ = V_cage_;
    SparseMatrix<double> S(num_vertices_cage_,num_vertices_cage_);// 稀疏矩阵表征顶点和顶点之间的相互联系强度？
    S.setIdentity();
    BiCGSTAB< SparseMatrix<double> > solverImpl;// 只是一个求解稀疏线性代数方程组的算子而已
    if (implicit_smoothing){
        S = S - lambda_smooth_ * laplace_beltrami_matrix_;
        solverImpl.compute(S);
    }
    else 
        S = S + lambda_smooth_ * laplace_beltrami_matrix_;

    Vector3d s, dir, dir_norm;
    vector<igl::Hit> hits;
    const double delta_T = 0.1;
    for (int i = 0; i < num_smoothing_iterations_; i++){// 光滑的迭代次数
        if (implicit_smoothing) {
            V_cage_smooth_ =  solverImpl.solve(V_cage_prev);
        } else {
            V_cage_smooth_ = S * V_cage_prev;
            // V_cage_smooth_ = (S * V_cage_prev.transpose()).transpose();
        }

        // Test for intersection with the mesh
        for (int j = 0; j < num_vertices_cage_; j++) {
            hits.clear();
            s = V_cage_prev.row(j);// 光滑前的顶点坐标
            dir = V_cage_smooth_.row(j) - V_cage_prev.row(j);// 同一个顶点 光滑后 - 光滑前 的坐标 (移动的delta)
            dir_norm = dir.normalized();
            if (ray_mesh_intersect(s, dir_norm, V_pca_basis_, F_mesh_, hits)){
                if (hits[0].t <= dir.norm()) {
                    V_cage_smooth_.row(j) = Vector3d(V_cage_prev.row(j)) + delta_T * hits[0].t * dir_norm;
                }
            }
        }
        V_cage_prev = V_cage_smooth_;
    }
}

// 返回cage的顶点在世界空间中的坐标
MatrixXd CageGenerator::GetCage() {
    MatrixXd result = (pca_basis_matrix_ * V_cage_.transpose()).transpose();
    for (int i = 0 ; i < result.rows(); i++){
        result.row(i) = barycenter_ + Vector3d(result.row(i));// 要加上barycenter_(整个mesh的中心)做平移
    }
    return result;
}

// 返回光滑后的cage的顶点在世界空间中的坐标
MatrixXd CageGenerator::GetSmoothCage() {
    MatrixXd result = (pca_basis_matrix_ * V_cage_smooth_.transpose()).transpose();
    for (int i = 0 ; i < result.rows(); i++){
        result.row(i) = barycenter_ + Vector3d(result.row(i));// 要加上barycenter_(整个mesh的中心)做平移
    }
    return result;
}