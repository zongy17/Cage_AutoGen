#include <iostream>

#include <igl/file_exists.h>
#include <time.h>
#include "CageGenerator.h"
#include "glRender.h"

using namespace Eigen;

float sparseness_cage = 0.5; // For automatic cage generation : CHANGE THIS for a sparser or denser cage : the larger the parameter the denser the cage

int main(int agrc, char* argv[]) 
{
    // LOAD MESHES
    MatrixXd V_mesh;
    MatrixXi F_mesh;

    std::string mesh_file_name = std::string(argv[1]);
    if (mesh_file_name.find("obj") != string::npos) 
        igl::readOBJ(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("ply") != string::npos) 
        igl::readPLY(mesh_file_name,V_mesh,F_mesh);
    else if (mesh_file_name.find("off") != string::npos) 
        igl::readOFF(mesh_file_name,V_mesh,F_mesh);
    else {
        cout << "Mesh file not recognized " << endl;
        return 1;
    }

    float lambda_smooth_implicit = .6;// 光滑系数
    int num_iterations_smoothing = 2;// 光滑迭代次数
    CageGenerator cage_generator(V_mesh, F_mesh, num_iterations_smoothing, lambda_smooth_implicit, sparseness_cage);

    MatrixXd bb_vertices;// 最粗的 bounding box 的顶点
    MatrixXi bb_faces;// 最粗的 bounding box 的面

    clock_t start = clock();
    cage_generator.ComputeCage(bb_vertices, bb_faces);// 自动生成cage
    std::cout << "Time to compute cage : " << (clock() - start) / (double) CLOCKS_PER_SEC << " sec" << std::endl << std::endl;
        
    MatrixXd V_cage_automatic = cage_generator.GetCage();
    MatrixXd V_cage_automatic_smooth = cage_generator.GetSmoothCage();
    MatrixXi F_cage_automatic = cage_generator.GetCageFaces();
    // 写出obj文件
    igl::writeOBJ("origin.obj", V_mesh, F_mesh);
    igl::writeOBJ("course.obj", V_cage_automatic, F_cage_automatic);
    igl::writeOBJ("smooth.obj", V_cage_automatic_smooth, F_cage_automatic);

    // renderLoop("course.obj", "smooth.obj");
    renderLoop("origin.obj", "course.obj", "smooth.obj");

    return 0;
}


