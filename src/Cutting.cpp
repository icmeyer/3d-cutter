#include <iostream>
#include <vector>
#include <array>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "Cutting.h"
#include <algorithm>
#include <map>
#include <fstream>
#include <sstream>


using namespace std;

#define d 3
#define T float
#define TV array<T, d>
#define T2 array<T, 2>
#define T4 array<T, 4>
#define Idp1 array<int, d+1>
#define Id array<int, d>
#define I3 array<int, 3>
#define I2 array<int, 2>
#define I4 array<int, 4>

#ifdef far
#undef far
#endif

#ifdef near
#undef near
#endif

TetMesh<T> tetMesh;
TriMesh<T> triMesh;

void simpleTet1Old() {
    tetMesh.nodes_.clear();
    tetMesh.mesh_.clear();
    tetMesh.nodes_.push_back(TV{0,0,0});
    tetMesh.nodes_.push_back(TV{1,0,0});
    tetMesh.nodes_.push_back(TV{0,1,0});
    tetMesh.nodes_.push_back(TV{0,0,1});
    tetMesh.nodes_.push_back(TV{1,1,1});
    tetMesh.mesh_.push_back(Idp1{0,1,2,3});
    tetMesh.mesh_.push_back(Idp1{4,2,1,3});
//    array<T,4> w;
//    Cutter3D<T>::computeIntersection(array<TV,4>{tetMesh.nodes_[0],tetMesh.nodes_[1],tetMesh.nodes_[2],tetMesh.nodes_[3]}, array<TV,1>{TV{0.1,0.1,0.1}}, w);
//    cout << w[0] << ", " << w[1] << ", " << w[2] << ", " << w[3] << endl;
    tetMesh.initializeSurfaceMesh();
    tetMesh.computeConnectedComponents();
    
    triMesh.nodes_.push_back(TV{0.3f,-0.1f,1.1f});
    triMesh.nodes_.push_back(TV{0.3f,-0.1f,-1});
    triMesh.nodes_.push_back(TV{0.3f,1.1f,1.1f});
    triMesh.nodes_.push_back(TV{0.3f,1.1f,-1});
    triMesh.mesh_.push_back(I3{0,1,2});
    triMesh.mesh_.push_back(I3{2,3,1});
    
    Cutter3D<T> cutter(tetMesh, triMesh);
    tetMesh = cutter.run();
    print<float, 3>(tetMesh.nodes_);
    print<int, 4>(tetMesh.mesh_);
    for (int &v: tetMesh.connectedComponents_){
        cout << v << endl;
    }
}

void simpleTet1() {
    vector<array<float, 3>> volumeNodes = {{0,0,0},
                                           {1,0,0},
                                           {0,1,0},
                                           {0,0,1},
                                           {1,1,1}};
    vector<array<int, 4>> volumeVoxels = {{0,1,2,3},
                                          {4,2,1,3}};

    tetMesh.initializeSurfaceMesh();
    tetMesh.computeConnectedComponents();
    vector<array<float, 3>> surfNodes = {{0.3, -0.1, 1.1},
                                         {0.3, -0.1, -1},
                                         {0.3, 1.1, 1.1},
                                         {0.3, 1.1, -1}};
    vector<array<int, 3>> surfFaces = {{0,1,2},
                                       {2,3,1}};
                                         
    Cutter3D<T> cutter(volumeNodes, volumeVoxels, surfNodes, surfFaces);
    cutter.runPlainOut();
    print<float, 3>(*cutter.outNodes);
    print<int, 4>(*cutter.outVoxels);
    for (int &v: *cutter.outIDs){
        cout << v << endl;
    }
      
}

void loadTriMesh(const string& filename) {
    ifstream fs(filename);
    triMesh.nodes_.clear();
    triMesh.mesh_.clear();
    string line;
    int l = 0;
    while(getline(fs, line)) {
        if (line.find("end_header") != std::string::npos) {
            while(getline(fs, line)) {
                if (line.substr(0,2) != "3 ") {
                    stringstream ss(line);
                    T t1,t2,t3;
                    ss >> t1 >> t2 >> t3;
                    triMesh.nodes_.push_back(divide<T,3>(TV{t1,t2,t3},13));
                } else {
                    stringstream ss(line);
                    int t1,t2,t3;
                    ss >> t1 >> t1 >> t2 >> t3;
                    triMesh.mesh_.push_back(I3{t1,t2,t3});
                }
            }
            break;
        }
    }
    cout << triMesh.nodes_.size() << " " << triMesh.mesh_.size() << endl;
}

void meshCutGrid() {
    tetMesh.nodes_.clear();
    tetMesh.mesh_.clear();
    int width = 23;
    int height = 23;
    int depth = 23;
    T low = -0.5;
    T high = 0.5;
    T left = -0.5;
    T right = 0.5;
    T far = 0.6;
    T near = -0.5;
    T lw = (right-left)/width;
    T lh = (high-low)/height;
    T ld = (far-near)/depth;
    int offset = (width+1)*(depth+1);

    for (int i = 0; i < height+1; i++) {
        //cout << "height: " << low+i*lh << endl;
        for (int j = 0; j < width+1; j++) {
            for (int k = 0; k < depth+1; k++) {
                tetMesh.nodes_.push_back(TV{left+j*lw, low+i*lh,near+k*ld});
            }
        }
    }
    
    int tetIndex = 0;
    for (int k = 0; k < height; k++) {
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < depth; j++) {
                tetMesh.mesh_.push_back(I4{tetIndex+offset+1, tetIndex+1, tetIndex, tetIndex+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex, tetIndex+offset, tetIndex+offset+1, tetIndex+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+offset, tetIndex+offset+1, tetIndex+offset+depth+1});
                tetMesh.mesh_.push_back(I4{tetIndex+1, tetIndex+depth+1, tetIndex+offset+1, tetIndex+depth+2});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+depth+2, tetIndex+offset+depth+2, tetIndex+offset+1});
                tetMesh.mesh_.push_back(I4{tetIndex+depth+1, tetIndex+depth+2+offset, tetIndex+offset+depth+1, tetIndex+offset+1});
                tetIndex++;
            }
            tetIndex++;
        }
        tetIndex += (depth+1);
    }

    tetMesh.initializeSurfaceMesh();
    tetMesh.computeConnectedComponents();
//    triMesh.nodes_.push_back(TV{0.3,-0.1,1.1});
//    triMesh.nodes_.push_back(TV{0.3,-0.1,-1});
//    triMesh.nodes_.push_back(TV{0.3,1.1,1.1});
//    triMesh.nodes_.push_back(TV{0.3,1.1,-1});
//    triMesh.mesh_.push_back(I3{0,1,2});
//    triMesh.mesh_.push_back(I3{2,3,1});

//    loadTriMesh("/Users/yutingwang/Downloads/beethoven.ply");
//    Cutter3D<T> cutter;
//    tetMesh = cutter.run(tetMesh, triMesh);
}


int main(int argc, char **argv) {
    if (argc==1){
        meshCutGrid();
        simpleTet1();
    }
    else if (argc==5){
        std::string tetNodes = argv[1];
        std::string tetTets = argv[2];
        std::string surfNodes = argv[3];
        std::string surfTris = argv[4];
        std::cout << "Running 3dcutter with: \n" 
                  << "Volume nodes:      " << tetNodes << "\n"
                  << "Volume tetrahedra: " << tetTets << "\n"
                  << "Surf nodes:        " << surfNodes << "\n"
                  << "Surf triangles:    " << surfTris << "\n";
    }
    return 0;
}
