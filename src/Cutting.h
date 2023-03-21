//
//  Copyright (c) 2018 Yuting Wang. All rights reserved.
//
#ifndef Cutting_h
#define Cutting_h

#include <vector>
#include <array>
#include <limits>
#include <algorithm>
#include <iostream>
#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <stack>
#include <map>
#include <numeric>
#include <math>


// Datastructures start

template<typename T, int d>
void print(const std::array<T, d>& v) {
    std::cout << "{";
    for (int i = 0; i < d; ++i) {
        std::cout << v[i] << ", ";
    }
    std::cout << "}\n";
}

template<typename T, int d>
void print(const std::vector<std::array<T, d>>& v) {
    std::cout << "{\n";
    for (const auto& a: v) {
        print<T,d>(a);
    }
    std::cout << "}\n";
}

template<typename T, int d1, int d2>
void print(const std::array<std::array<T, d1>, d2>& v) {
    std::cout << "{\n";
    for (const auto& a: v) {
        print<T,d1>(a);
    }
    std::cout << "}\n";
}

template<typename T>
T sorted(const T& v) {
    auto sv = v;
    std::sort(sv.begin(), sv.end());
    return sv;
}

class UnionFind {
    int *id, cnt, *sz;
public:
    // Create an empty union find data structure with N isolated sets.
    UnionFind(int N)   {
        cnt = N;
        id = new int[N];
        sz = new int[N];
        for(int i=0; i<N; i++)	{
            id[i] = i;
            sz[i] = 1;
        }
    }
    ~UnionFind()	{
        delete [] id;
        delete [] sz;
    }
    // Return the id of component corresponding to object p.
    int find(int p)	{
        int root = p;
        while (root != id[root])
            root = id[root];
        while (p != root) {
            int newp = id[p];
            id[p] = root;
            p = newp;
        }
        return root;
    }
    // Replace sets containing x and y with their union.
    void merge(int x, int y)	{
        int i = find(x);
        int j = find(y);
        if (i == j) return;
        
        // make smaller root point to larger one
        if   (sz[i] < sz[j])	{
            id[i] = j;
            sz[j] += sz[i];
        } else	{
            id[j] = i;
            sz[i] += sz[j];
        }
        cnt--;
    }
    // Are objects x and y in the same set?
    bool connected(int x, int y)    {
        return find(x) == find(y);
    }
    // Return the number of disjoint sets.
    int count() {
        return cnt;
    }
};

template<typename T, int d>
std::array<T,4> toI4(const std::array<T,d>& Id, T fill = -1) {
    std::array<T,4> a;
    a.fill(fill);
    for (int i = 0; i < std::min(d, 4); ++i) {
        a[i] = Id[i];
    }
    return a;
}
template<typename T, int d>
std::array<T, d> add(const std::array<T, d>& a1, const std::array<T, d>& a2) {
    std::array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a1[j] + a2[j];
    }
    return s;
}

template<typename T, int d>
std::array<T, d> substract(const std::array<T, d>& a1, const std::array<T, d>& a2) {
    std::array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a1[j] - a2[j];
    }
    return s;
}

template<typename T, int d>
std::array<T, d> divide(const std::array<T, d>& a, T t) {
    std::array<T, d> s;
    for (int j = 0; j < d; ++j) {
        s[j] = a[j] / t;
    }
    return s;
}

template<typename T, int d, int n>
std::array<T, d> center(const std::array<std::array<T, d>, n>& nodes, const std::array<T, n> weights) {
    std::array<T, d> c;
    for (size_t i = 0; i < d; ++i) {
        c[i] = 0;
        for (size_t j = 0; j < n; ++j) {
            c[i] += nodes[j][i] * weights[j];
        }
    }
    return c;
}

template<typename T, int d>
std::array<T, d> elementCenter(const std::vector<std::array<T,d>>& nodes, const std::array<int, 4>& element, const std::array<T, 4>& weights) {
    std::array<T, d> c;
    c.fill(0);
    for (size_t i = 0; i < 4; ++i) {
        if (element[i] >= 0) {
            const auto& node = nodes[element[i]];
            auto w = weights[i];
            for (int j = 0; j < d; ++j) {
                c[j] += (node[j] * w);
            }
        } else {
            break;
        }
    }
    return c;
}

template<typename T, int d>
T dot(const std::array<T, d>& a1, const std::array<T, d>& a2) {
    T s = 0;
    for (int j = 0; j < d; ++j) {
        s += a1[j] * a2[j];
    }
    return s;
}

template<typename T>
T cross(const std::array<T, 2>& a1, const std::array<T, 2>& a2) {
    return a1[0] * a2[1] - a1[1] * a2[0];
}

template<typename T>
std::array<T, 3> cross(const std::array<T, 3>& a1, const std::array<T, 3>& a2) {
    std::array<T, 3> s;
    for (int j = 0; j < 3; ++j) {
        s[j] = a1[(j+1)%3] * a2[(j+2)%3] - a2[(j+1)%3] * a1[(j+2)%3];
    }
    return s;
}

template<typename T>
T volume(const std::array<T,3>& node1, const std::array<T,3>& node2, const std::array<T,3>& node3, const std::array<T,3>& node4) {
    return dot<T,3>(cross<T>(substract<T,3>(node2,node1), substract<T,3>(node3,node1)), substract<T,3>(node4,node1));
}

template<typename T, int d>
T norm(const std::array<T, d>& a1) {
    return sqrt(dot(a1, a1));
}

template<typename T, int d>
T pointEdgeWeight(const std::array<T, d>& e1, const std::array<T, d>& e2, const std::array<T, d>& p) {
    std::array<T, d> v1 = substract<T, d>(p, e1);
    std::array<T, d> v2 = substract<T, d>(e2, e1);
    return dot<T, d>(v1, v2) / dot<T, d>(v2, v2);
}

template<typename T, int d>
bool pointOnEdge(const std::array<T, d>& e1, const std::array<T, d>& e2, const std::array<T, d>& p, T& w) {
    std::array<T, d> v1 = substract<T, d>(p, e1);
    std::array<T, d> v2 = substract<T, d>(e1, e2);
    if (cross<T>(v1, v2) == 0) {
        w = pointEdgeWeight<T, d>(e1, e2, p);
        return true;
    }
    return false;
}

template<typename T>
bool edgeEdgeIntersect(const std::array<T, 2>& p1, const std::array<T, 2>& p2, const std::array<T, 2>& q1, const std::array<T, 2>& q2, T& w) {
    std::array<T, 2> e1 = substract<T, 2>(p2, p1);
    T a1 = cross<T>(substract<T, 2>(q1, p1), e1);
    T a2 = cross<T>(substract<T, 2>(q2, p1), e1);
    if ((a1 < 0 && a2 > 0) || (a1 > 0 && a2 < 0)) {
        std::array<T, 2> e2 = substract<T, 2>(q2, q1);
        T a3 = cross<T>(substract<T, 2>(p1, q1), e2);
        T a4 = cross<T>(substract<T, 2>(p2, q1), e2);
        if ((a3 < 0 && a4 > 0) || (a3 > 0 && a4 < 0)) {
            w = a3 / (a3 - a4);
            return true;
        }
    }
    return false;
}

template<typename T>
bool pointInTriangle(const std::array<std::array<T, 2>, 3>& triangle, const std::array<T, 2>& point, std::array<T, 3>& w) {
    std::array<T, 3> areas;
    for (int i = 0; i < 3; ++i) {
        areas[(i+2)%3] = cross<T>(substract<T, 2>(point, triangle[i]), substract<T, 2>(triangle[(i+1)%3], triangle[i]));
    }
    if ((areas[0] < 0 && areas[1] < 0 && areas[2] < 0) || (areas[0] > 0 && areas[1] > 0 && areas[2] > 0)) {
        w = divide<T, 3>(areas, areas[0]+areas[1]+areas[2]);
        return true;
    }
    return false;
}

template<typename T, int d, int d2>
std::array<T, d> elementCenter(const std::vector<std::array<T, d>>& vertices, const std::array<int, d2>& element) {
    std::array<T, d> center;
    for (auto i : element) {
        for (int j = 0; j < d; ++j) {
            center[j] += vertices[i][j];
        }
    }
    for (int i = 0; i < d; ++i) {
        center[i] /= (T)element.size();
    }
    return center;
}

template<typename T, int d>
struct Box {
    std::array<T, d> lowerLeft_, upperRight_;
    
    Box() {}
    
    Box(const Box& b) {
        for (int i = 0; i < d; ++i) {
            lowerLeft_[i] = b.lowerLeft_[i];
            upperRight_[i] = b.upperRight_[i];
        }
    }
    
    Box(const Box& b1, const Box& b2) {
        for (int i = 0; i < d; ++i) {
            lowerLeft_[i] = std::min(b1.lowerLeft_[i], b2.lowerLeft_[i]);
            upperRight_[i] = std::max(b1.upperRight_[i], b2.upperRight_[i]);
        }
    }
    
    bool intersects(const Box& b) {
        for (int i = 0; i < d; ++i) {
            if (lowerLeft_[i] > b.upperRight_[i] || upperRight_[i] < b.lowerLeft_[i]) {
                return false;
            }
        }
        return true;
    }
};

template<typename T, int d, int d1>
Box<T,d> buildBox(const std::vector<std::array<T, d>>& vertices, const std::array<int, d1>& element) {
    Box<T,d> b;
//    print<int,d1>(element);
    for (int i = 0; i < d; ++i) {
        b.lowerLeft_[i] = std::numeric_limits<T>::max();
        b.upperRight_[i] = std::numeric_limits<T>::lowest();
    }
//    print<T,d>(b.lowerLeft_);
//    print<T,d>(b.upperRight_);
    for (size_t i = 0; i < d1; ++i) {
//        print<T,d>(vertices[element[i]]);
        for (int j = 0; j < d; ++j) {
            b.lowerLeft_[j] = std::min(b.lowerLeft_[j], vertices[element[i]][j]);
            b.upperRight_[j] = std::max(b.upperRight_[j], vertices[element[i]][j]);
        }
    }
//    print<T,d>(b.lowerLeft_);
//    print<T,d>(b.upperRight_);
    return b;
}

template<typename T, int d>
struct BoxNode {
    int n_;
    BoxNode *left_, *right_;
    Box<T, d> box_;
    
    BoxNode(int n, const std::vector<Box<T, d>>& boxes) : n_(n), left_(nullptr), right_(nullptr), box_(boxes[n]) {}
    
    BoxNode(BoxNode* left, BoxNode* right) : n_(-1), left_(left), right_(right), box_(left->box_, right->box_) {}
    
    ~BoxNode() {
        if (left_ != nullptr) {
            delete left_;
        }
        if (right_ != nullptr) {
            delete right_;
        }
    }
};

template<typename T, int d>
BoxNode<T, d>* buildBoxHierarchy(const std::vector<Box<T, d>>& boxes, const std::vector<std::array<T, d>>& centers, std::vector<size_t>& elementIndexes, int begin, int end, int level) {
    BoxNode<T, d>* root = nullptr;
    if (elementIndexes.size() == 0) {
        return nullptr;
    }
    if (begin == end) {
        root = new BoxNode<T, d>(elementIndexes[begin], boxes);
    } else {
        nth_element(elementIndexes.begin()+begin, elementIndexes.begin()+(begin+end)/2, elementIndexes.begin()+end, [&](std::size_t i, std::size_t j){ return centers[i][level%d] < centers[j][level%d]; });
        BoxNode<T, d> *left = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, begin, (begin+end)/2, ++level);
        BoxNode<T, d> *right = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, (begin+end)/2+1, end, level);
        root = new BoxNode<T, d>(left, right);
    }
    return root;
}

template<typename T, int d>
class BoxHierarchy {
    int n_; //number of boxes
    BoxNode<T, d>* root_;
public:
    BoxHierarchy(const std::vector<Box<T, d>>& boxes, const std::vector<std::array<T, d>>& centers) : n_(centers.size()) {
        if (boxes.size()) {
            std::vector<size_t> elementIndexes(boxes.size());
            iota(elementIndexes.begin(), elementIndexes.end(), 0);
            root_ = buildBoxHierarchy<T, d>(boxes, centers, elementIndexes, 0, boxes.size()-1, 0);
        }
    }
    
    void intersect(const BoxHierarchy<T, d>& bh, std::vector<std::vector<int>>& intersectingElements) const {
        intersectingElements.clear();
        intersectingElements.resize(n_);
        if (!root_ || !bh.root_) {
            return;
        }
        std::stack<std::pair<BoxNode<T, d>*, BoxNode<T, d>*>> s;
        s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(root_, bh.root_));
        while (s.size()) {
            std::pair<BoxNode<T, d>*, BoxNode<T, d>*> top = s.top();
            s.pop();
            if (top.first->box_.intersects(top.second->box_)) {
                if (top.first->n_ != -1 && top.second->n_ != -1) {
                    intersectingElements[top.first->n_].push_back(top.second->n_);
                } else if (top.second->n_ != -1) {
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second));
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second));
                } else if (top.first->n_ != -1) {
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first, top.second->left_));
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first, top.second->right_));
                } else {
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second->left_));
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->left_, top.second->right_));
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second->left_));
                    s.push(std::pair<BoxNode<T, d>*, BoxNode<T, d>*>(top.first->right_, top.second->right_));
                }
            }
        }
    }
    
    ~BoxHierarchy() {
        delete root_;
    }
};

template<typename T, int d1, int d2>
BoxHierarchy<T, d1> buildBoxHierarchy(const std::vector<std::array<T,d1>>& nodes, const std::vector<std::array<int, d2>>& elements) {
    std::vector<Box<T, d1>> boxes;
    std::vector<std::array<T, d1>> centers;
    for (const auto& e: elements) {
        boxes.push_back(buildBox<T,d1,d2>(nodes, e));
        centers.push_back(elementCenter<T,d1,d2>(nodes, e));
    }
    return BoxHierarchy<T, d1>(boxes, centers);
}

template<typename T>
class TriMesh {
    typedef std::array<int, 2> I2;
    typedef std::array<int, 3> I3;
    typedef std::array<int, 4> I4;
    typedef std::array<T, 3> TV;
    
public:
    std::vector<TV> nodes_;
    std::vector<I3> mesh_;
    
    void clear() {
        nodes_.clear();
        mesh_.clear();
    }
};

const std::array<std::array<int, 3>, 4> FaceIndexes = {
    std::array<int, 3>{0,1,2},
    std::array<int, 3>{0,2,3},
    std::array<int, 3>{1,2,3},
    std::array<int, 3>{0,3,1}
};

std::array<std::array<int,3>,4> tetFaces(const std::array<int,4>& tet) {
    std::array<std::array<int,3>,4> faces;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            faces[i][j] = tet[FaceIndexes[i][j]];
        }
    }
    return faces;
}

std::array<int,3> tetFace(const std::array<int,4>& tet, int i) {
    std::array<int,3> face;
    for (int j = 0; j < 3; ++j) {
        face[j] = tet[FaceIndexes[i][j]];
    }
    return face;
}

const std::array<std::array<int, 2>, 6> EdgeIndexes = {std::array<int,2>{0,1}, std::array<int,2>{0,2}, std::array<int,2>{0,3}, std::array<int,2>{1,2}, std::array<int,2>{1,3}, std::array<int,2>{2,3}};

std::array<std::array<int,2>,6> tetEdges(const std::array<int,4>& tet) {
    std::array<std::array<int,2>,6> faces;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 2; ++j) {
            faces[i][j] = tet[EdgeIndexes[i][j]];
        }
    }
    return faces;
}

std::array<std::array<int,2>,3> faceEdges(const std::array<int,3>& face) {
    std::array<std::array<int,2>,3> edges;
    for (int i = 0; i < 3; ++i) {
        edges[i] = std::array<int,2>{face[i], face[(i+1)%3]};
    }
    return edges;
}

template<typename T, int d, int d1>
std::array<std::array<T,d>,d1> elementNodes(const std::vector<std::array<T,d>>& nodes, const std::array<int,d1>& element) {
    std::array<std::array<T,d>,d1> ps;
    for (int i = 0; i < d1; ++i) {
        ps[i] = nodes[element[i]];
    }
    return ps;
}


template<typename T>
class TetMesh {
    typedef std::array<int, 2> I2;
    typedef std::array<int, 3> I3;
    typedef std::array<int, 4> I4;
    typedef std::array<T, 3> TV;

public:
    std::vector<TV> nodes_;
    std::vector<I4> mesh_;
    std::vector<I3> surfaceMesh_;
    std::vector<int> connectedComponents_; //connected component id of each element

    TetMesh() {}

    TetMesh(std::vector<TV>&& nodes, std::vector<I4>&& mesh): nodes_(nodes), mesh_(mesh) {
        initializeSurfaceMesh();
        computeConnectedComponents();
    }
    void initializeSurfaceMesh() {
        surfaceMesh_.clear();
        std::map<I3, I3> surfaceElements; //sorted to unsorted elements
        for (const auto& tet: mesh_) {
            for (const auto& fi: FaceIndexes) {
                auto face = I3{tet[fi[0]], tet[fi[1]], tet[fi[2]]};
                auto sortedFace = face;
                std::sort(sortedFace.begin(), sortedFace.end());
                if (surfaceElements.count(sortedFace)) {
                    surfaceElements.erase(sortedFace);
                } else {
                    surfaceElements[sortedFace] = face;
                }
            }
        }
        for (const auto& e: surfaceElements) {
            surfaceMesh_.push_back(e.second);
        }
    }
    
    void computeConnectedComponents() {
        connectedComponents_.resize(mesh_.size());
        UnionFind nodeClasses(nodes_.size());
        for (int i = 0; i < mesh_.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                nodeClasses.merge(mesh_[i][j], mesh_[i][j+1]);
            }
        }
        std::map<int, int> nodeClassToCC;
        int c = 1;
        for (int i = 0; i < mesh_.size(); ++i) {
            if (!nodeClassToCC.count(nodeClasses.find(mesh_[i][0]))) {
                connectedComponents_[i] = c;
                nodeClassToCC[nodeClasses.find(mesh_[i][0])] = c;
                ++c;
            } else {
                connectedComponents_[i] = nodeClassToCC[nodeClasses.find(mesh_[i][0])];
            }
        }
        std::cout << "found " << c-1 << " connected components\n";
    }
};

// Datastructures end

template<typename T>
class Cutter3D {
    typedef std::array<int, 1> I1;
    typedef std::array<int, 2> I2;
    typedef std::array<int, 3> I3;
    typedef std::array<int, 4> I4;
    typedef std::array<int, 5> I5;
    typedef std::array<T, 2> T2;
    typedef std::array<T, 3> T3;
    typedef std::array<T, 4> T4;
    typedef std::map<I4, T4> Intersections;
    typedef std::map<I4, std::vector<int>> TetBoundary2TetIds;

    struct CutElement {
        int parentElementIndex;
        std::array<bool, 4> subElements; // in the same order as the tet nodes
        
        CutElement(int i, bool fill = true): parentElementIndex(i) {
            subElements.fill(fill);
        }

        int numPieces() const {
            return (int)subElements[0] + (int)subElements[1] + (int)subElements[2] + (int)subElements[3];
        }
    };
    
    static bool computeIntersection(const std::array<T3,2>& nodes1, const std::array<T3,3>& nodes2, std::array<T,2>& w1, std::array<T, 3>& w2) {
        T v1 = volume<T>(nodes1[0], nodes2[0], nodes2[1], nodes2[2]);
        T v2 = volume<T>(nodes1[1], nodes2[0], nodes2[1], nodes2[2]);
        T v3 = volume<T>(nodes1[0], nodes1[1], nodes2[0], nodes2[1]);
        T v4 = volume<T>(nodes1[0], nodes1[1], nodes2[1], nodes2[2]);
        T v5 = volume<T>(nodes1[0], nodes1[1], nodes2[2], nodes2[0]);
        if (v1*v2<0 && (v3>0)==(v4>0) && (v4>0)==(v5>0)) {
            w1[0] = fabs(v2) / (fabs(v1) + fabs(v2));
            w1[1] = 1 - w1[0];
            T v = fabs(v3) + fabs(v4) + fabs(v5);
            w2[0] = fabs(v4) / v;
            w2[1] = fabs(v5) / v;
            w2[2] = 1 - w2[0] - w2[1];
            //cout << w1[0] << std::endl;
            return true;
        } else {
            return false;
        }
    }

    static bool computeIntersection(const std::array<T3,2>& nodes1, const std::array<T3,3>& nodes2, std::array<T,2>& w) {
        std::array<T,3> w1;
        return computeIntersection(nodes1, nodes2, w, w1);
    }

    static bool computeIntersection(const std::array<T3,3>& nodes1, const std::array<T3,2>& nodes2, std::array<T,3>& w) {
        std::array<T,2> w1;
        return computeIntersection(nodes2, nodes1, w1, w);
    }

    static bool computeIntersection(const std::array<T3,4>& nodes1, const std::array<T3,1>& nodes2, std::array<T,4>& w) {
        T v1 = volume<T>(nodes1[0], nodes1[1], nodes1[2], nodes2[0]);
        T v2 = volume<T>(nodes1[0], nodes1[2], nodes1[3], nodes2[0]);
        T v3 = volume<T>(nodes1[0], nodes1[3], nodes1[1], nodes2[0]);
        T v4 = volume<T>(nodes2[0], nodes1[1], nodes1[2], nodes1[3]);
        if (v1 == 0 || v2 == 0 || v3 == 0 || v4 == 0) {
            std::cout << "point tet degenerate case" << std::endl;
        }
        // cout << v1 << ", " << v2 << ", " << v3 << ", " << v4 << std::endl;
        if ((v1>0) == (v2>0) && (v2>0) == (v3>0) && (v3>0) == (v4>0)) {
            T v = fabs(v1) + fabs(v2) + fabs(v3) + fabs(v4);
            w[0] = fabs(v4) / v;
            w[1] = fabs(v2) / v;
            w[2] = fabs(v3) / v;
            w[3] = 1 - w[0] - w[1] - w[2];
            return true;
        } else {
            return false;
        }
        return false;
    }

    template<int d1, int d2>
    static void computeIntersections(const std::vector<T3>& nodes1, const std::vector<T3>& nodes2, const std::vector<std::array<int,d1>>& e1, const std::vector<std::array<int,d2>>& e2, const BoxHierarchy<T,3>& b1, const BoxHierarchy<T,3>& b2, std::map<I4, T4>& intersections) {
        std::vector<std::vector<int>> intersectingBoxes; // intersecting boxes
        b1.intersect(b2, intersectingBoxes);
        for (size_t i = 0; i < intersectingBoxes.size(); ++i) {
            //cout << "e1 " << i << std::endl;
            //print<int,d1>(e1[i]);
            for (auto j : intersectingBoxes[i]) {
                //cout << j << ", " << std::endl;
                //print<int,d2>(e2[j]);
                auto tetNodes = elementNodes<T,3,d1>(nodes1, e1[i]);
                auto triNodes = elementNodes<T,3,d2>(nodes2, e2[j]);
                std::array<T,d1> w;
                if (computeIntersection(tetNodes, triNodes, w)) {
                    intersections[toI4<int,d1>(e1[i])] = toI4<T,d1>(w,0);
                }
            }
        }
    }

    static Intersections computeIntersections(const TetMesh<T>& tetMesh, const TriMesh<T>& triMesh, TetBoundary2TetIds& tetBoundary2TetIds) {
        std::map<I4, T4> intersections;
        
        // build box hierarchies for tetMesh
        std::set<I3> tetMeshFaces;
        std::set<I2> tetMeshEdges;
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            auto tet = tetMesh.mesh_[i];
            std::sort(tet.begin(), tet.end());
            tetBoundary2TetIds[tet].push_back(i);
            auto faces = tetFaces(tet);
            for (auto& face: faces) {
                std::sort(face.begin(), face.end());
                tetBoundary2TetIds[toI4<int,3>(face)].push_back(i);
                tetMeshFaces.insert(face);
            }
            auto edges = tetEdges(tet);
            for (auto& edge: edges) {
                std::sort(edge.begin(), edge.end());
                tetBoundary2TetIds[toI4<int,2>(edge)].push_back(i);
                tetMeshEdges.insert(edge);
            }
        }
        std::vector<I1> tetMeshNodeVec;
        for (int i = 0; i < tetMesh.nodes_.size(); ++i) {
            tetMeshNodeVec.push_back(I1{i});
        }
        std::vector<I3> tetMeshFaceVec(tetMeshFaces.begin(), tetMeshFaces.end());
        std::vector<I2> tetMeshEdgeVec(tetMeshEdges.begin(), tetMeshEdges.end());
        std::cout << "buliding tet mesh hierarchy" << std::endl;
        auto tetMeshHierarchy = buildBoxHierarchy<T,3,4>(tetMesh.nodes_, tetMesh.mesh_);
        auto tetMeshFaceHierarchy = buildBoxHierarchy<T,3,3>(tetMesh.nodes_, tetMeshFaceVec);
        auto tetMeshEdgeHierarchy = buildBoxHierarchy<T,3,2>(tetMesh.nodes_, tetMeshEdgeVec);
        auto tetMeshNodeHierarchy = buildBoxHierarchy<T,3,1>(tetMesh.nodes_, tetMeshNodeVec);
        std::cout << "tet mesh hierarchy built" << std::endl;

        // box hierarchy for triMesh
        std::set<I2> triMeshEdges;
        for (const auto& tri: triMesh.mesh_) {
            auto edges = faceEdges(tri);
            for (auto& edge: edges) {
                std::sort(edge.begin(), edge.end());
                triMeshEdges.insert(edge);
            }
        }
        std::vector<I2> triMeshEdgeVec(triMeshEdges.begin(), triMeshEdges.end());
        std::vector<I1> triMeshNodeVec;
        for (int i = 0; i < triMesh.nodes_.size(); ++i) {
            triMeshNodeVec.push_back(I1{i});
        }
//        std::cout << "trimeshhierarchy\n";
//        print<T,3>(triMesh.nodes_);
//        print<int,3>(triMesh.mesh_);
        auto triMeshHierarchy = buildBoxHierarchy<T,3,3>(triMesh.nodes_, triMesh.mesh_);
//        std::cout << "trimeshhierarchy\n";
        auto triMeshEdgeHierarchy = buildBoxHierarchy<T,3,2>(triMesh.nodes_, triMeshEdgeVec);
        auto triMeshNodeHierarchy = buildBoxHierarchy<T,3,1>(triMesh.nodes_, triMeshNodeVec);
        std::cout << "tri mesh hierarchy built" << std::endl;

        // compute intersections
        // v-v
        // v-e
        // v-f
        // e-v
        // e-e
        // e-f
        computeIntersections<2,3>(tetMesh.nodes_, triMesh.nodes_, tetMeshEdgeVec, triMesh.mesh_, tetMeshEdgeHierarchy, triMeshHierarchy, intersections);
        // f-v
        // f-e
        computeIntersections<3,2>(tetMesh.nodes_, triMesh.nodes_, tetMeshFaceVec, triMeshEdgeVec, tetMeshFaceHierarchy, triMeshEdgeHierarchy, intersections);
        // t-v
        computeIntersections<4,1>(tetMesh.nodes_, triMesh.nodes_, tetMesh.mesh_, triMeshNodeVec, tetMeshHierarchy, triMeshNodeHierarchy, intersections);

        return intersections;
    }

    static std::vector<CutElement> split(const TetMesh<T>& tetMesh, const Intersections& intersections, TetBoundary2TetIds& tetBoundary2TetIds, std::set<int>& cutTets) {
        cutTets.clear();
        for (const auto& t: tetBoundary2TetIds) {
            if (intersections.count(t.first)) {
                for (auto i: t.second) {
                    cutTets.insert(i);
                }
            }
        }
        std::cout << cutTets.size() << " tets cut\n";
        std::vector<CutElement> v;
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            if (cutTets.count(i)) {
                std::array<bool,4> added;
                added.fill(false);
                auto tet = tetMesh.mesh_[i];
                for (int j = 0; j < 4; ++j) {
                    if (!added[j]) {
                        // find all connected pieces
                        CutElement ce(i, false);
                        std::stack<int> s;
                        s.push(j);
                        while(s.size()) {
                            auto top = s.top();
                            ce.subElements[top] = true;
                            added[top] = true;
                            s.pop();
                            // add all the connected pieces that are not added yet
                            for (int k = 0; k < 4; ++k) {
                                if (!added[k]) {
                                    if (!intersections.count(toI4<int,2>(sorted(I2{tet[top],tet[k]})))) {
                                        s.push(k);
                                    }
                                }
                            }
                        }
                        v.push_back(ce);
                    }
                }
            }
        }
        return v;
    }

    void static newTet(int parentId, const I4& tet, const TetMesh<T>& tetMesh, std::vector<T3>& newNodes, std::vector<I4>& newMesh, std::map<int,int>& nodeMapping, UnionFind& uf) {
        I4 newTet;
        //cout << "parent id " << parentId << std::endl;
        for (int i = 0; i < 4; ++i) { // for each node
            int newId = uf.find(tet[i]);
            //cout << tet[i] << ", " << newId << std::endl;
            const auto& it = nodeMapping.find(newId);
            if (it != nodeMapping.end()) {
                newTet[i] = it->second;
            } else {
                newTet[i] = newNodes.size();
                nodeMapping[newId] = newNodes.size();
                newNodes.push_back(tetMesh.nodes_[tetMesh.mesh_[parentId][i]]);
            }
        }

        newMesh.push_back(newTet);
    }
    
    static void merge(const std::vector<CutElement>& cutElements, const TetMesh<T>& tetMesh, std::vector<T3>& newNodes, std::vector<I4>& newMesh, const Intersections& intersections) {
        newNodes.clear();
        newMesh.clear();
        UnionFind uf(tetMesh.nodes_.size() + 4 * cutElements.size());
        std::map<I5, int> faceNode2NewNode; // key = {face,materialNode,node}
        std::set<int> cutTets;
        int total = tetMesh.nodes_.size();
        for (const auto& ce: cutElements) { // need to do face-face merging even for tets that are touched by the cut but not split, so that if a neighbor splits they are all connected to it.
            cutTets.insert(ce.parentElementIndex);
            const auto& tet = tetMesh.mesh_[ce.parentElementIndex];
            for (int i = 0; i < 4; ++i) { // for each face
                auto face = tetFace(tet, i);
                std::sort(face.begin(), face.end());
                I5 key;
                for (int j = 0; j < 3; ++j) {
                    key[j] = face[j];
                }
                for (int j = 0; j < 3; ++j) { // for each node check for material
                    int fij = FaceIndexes[i][j];
                    if (ce.subElements[fij]) {
                        key[3] = tet[fij];
                        uf.merge(total+fij, key[3]);
                        for (int k = 0; k < 3; ++k) { // for each node, merge
                            int fik = FaceIndexes[i][k];
                            key[4] = tet[fik];
                            int newId = total+fik;
                            //print<int,5>(key);
                            const auto& it = faceNode2NewNode.find(key);
                            if (it != faceNode2NewNode.end()) {
                                //cout << "merging " << it->second << ", " << newId << std::endl;
                                uf.merge(it->second, newId);
                            } else {
                                faceNode2NewNode[key] = newId;
                            }
                        }
                    }
                }
            }
            total += 4;
        }
        total = tetMesh.nodes_.size();
        std::map<int,int> nodeMapping;
        for (const auto& ce: cutElements) {
            newTet(ce.parentElementIndex, I4{total, total+1, total+2, total+3}, tetMesh, newNodes, newMesh, nodeMapping, uf);
            total += 4;
        }
        for (int i = 0; i < tetMesh.mesh_.size(); ++i) {
            if (!cutTets.count(i)) {
                newTet(i, tetMesh.mesh_[i], tetMesh, newNodes, newMesh, nodeMapping, uf);
            }
        }

//        std::cout << "merged mesh \n";
//        print<T,3>(newNodes);
//        print<int,4>(newMesh);
    }
    
    static TetMesh<T> subdivide(const std::vector<CutElement>& cutElements, const TetMesh<T>& tetMesh, std::vector<T3>& newNodes, std::vector<I4>& newMesh, Intersections& intersections) {
        // add a new node inside the tet, connect with cuts on each face to subdivide the tet
        std::map<I4, int> newNodeMapping;
        for (int i = 0; i < cutElements.size(); ++i) {
            const auto& ce = cutElements[i];
            const auto& originalTet = tetMesh.mesh_[ce.parentElementIndex];
            const auto sortedOriginalTet = sorted(originalTet);
            const auto& tet = newMesh[i];
            
            // get all edge cuts and add them as new nodes
            const auto originalEdges = tetEdges(originalTet);
            const auto edges = tetEdges(tet);
            int cutEdges = 0;
            T4 averageEdgeWeight{0,0,0,0};
            std::map<int, T> originalNodeId2Weight;
            for (int k = 0; k < originalEdges.size(); ++k) {
                auto sortedOriginalEdge = toI4<int,2>(sorted(originalEdges[k]));
                auto sortedEdge = toI4<int,2>(sorted(edges[k]));
                const auto& it = intersections.find(sortedOriginalEdge);
                if (it != intersections.end()) {
                    ++cutEdges;
                    for (int j = 0; j < 2; ++j) {
                        originalNodeId2Weight[sortedOriginalEdge[j]] += it->second[j];
                    }
                    const auto& idIt = newNodeMapping.find(sortedEdge);
                    if (idIt == newNodeMapping.end()) {
                        newNodeMapping[sortedEdge] = newNodes.size();
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalEdge, it->second));
//                        std::cout << "edge node ";
//                        print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalEdge, it->second));
                    }
                }
            }
            for (int j = 0; j < 4; ++j) {
                averageEdgeWeight[j] = originalNodeId2Weight[sortedOriginalTet[j]];
            }
            //cout << "cutEdges " << cutEdges << std::endl;

            // face cuts
            const auto originalFaces = tetFaces(originalTet);
            const auto faces = tetFaces(tet);
            for (int k = 0; k < faces.size(); ++k) {
                auto sortedOriginalFace = toI4<int,3>(sorted(originalFaces[k]));
                auto sortedFace = toI4<int,3>(sorted(faces[k]));
                const auto& it = intersections.find(sortedOriginalFace);
                if (it != intersections.end()) { // face center already computed
                    const auto& idIt = newNodeMapping.find(sortedFace);
                    if (idIt == newNodeMapping.end()) {
                        newNodeMapping[sortedFace] = newNodes.size();
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, it->second));
                    }
//                    std::cout << "face center ";
//                    print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, it->second));
                } else { // use average of edge cuts if not
                    int numEdges = 0;
                    T4 faceWeights{0,0,0,0};
                    std::map<int, T> node2weight;
                    for (int j = 0; j < 3; ++j) {
                        auto sortedOriginalEdge = toI4<int,2>(sorted(std::array<int,2>{sortedOriginalFace[j], sortedOriginalFace[(j+1)%3]}));
                        const auto& edgeIt = intersections.find(sortedOriginalEdge);
                        if (edgeIt != intersections.end()) {
                            ++numEdges;
                            for (int e = 0; e < 2; ++e) {
                                node2weight[sortedOriginalEdge[e]] += edgeIt->second[e];
                            }
                        }
                    }
                    if (numEdges > 1) { // otherwise don't add new face center
                        newNodeMapping[sortedFace] = newNodes.size();
                        for (int j = 0; j < 3; ++j) {
                            faceWeights[j] = node2weight[sortedOriginalFace[j]] / numEdges;
                        }
//                        std::cout << "face weight ";
//                        print<T,4>(faceWeights);
//                        std::cout << "face center ";
//                        print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, faceWeights));
                        newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalFace, faceWeights));
                        intersections[sortedOriginalFace] = faceWeights;
                    }
                }
            }
            
            // tet center
            int tetCenterId = newNodes.size();
            const auto& tetCenterIt = intersections.find(sortedOriginalTet);
            if (tetCenterIt != intersections.end()) {
                newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, tetCenterIt->second));
//                std::cout << "tet center ";
//                print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, tetCenterIt->second));
            } else { // if doesn't exist, use average of edge cuts or the center
                if (ce.numPieces() == 4) {
                    averageEdgeWeight.fill(0.25);
                } else {
                    averageEdgeWeight = divide<T,4>(averageEdgeWeight, cutEdges);
//                    print<T,4>(averageEdgeWeight);
                }
                newNodes.push_back(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, averageEdgeWeight));
//                std::cout << "tet center ";
//                print<T,3>(elementCenter<T,3>(tetMesh.nodes_, sortedOriginalTet, averageEdgeWeight));
                intersections[sortedOriginalTet] = averageEdgeWeight;
            }

            // add elements that are created by the new nodes added above
            std::vector<I4> newTets;
            for (int f = 0; f < faces.size(); ++f) {
                const auto& face = faces[f];
                const auto sortedFace = toI4<int,3>(sorted(face));
                const auto& newFaceCenterIt = newNodeMapping.find(sortedFace);
                if (newFaceCenterIt != newNodeMapping.end()) {
                    for (int j = 0; j < 3; ++j) {
                        auto sortedEdge = toI4<int,2>(sorted(std::array<int,2>{face[j], face[(j+1)%3]}));
                        const auto& newEdgeCenterIt = newNodeMapping.find(sortedEdge);
                        if (newEdgeCenterIt != newNodeMapping.end()) {
                            if (ce.subElements[FaceIndexes[f][j]]) {
                                newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, face[j], newEdgeCenterIt->second});
                            }
                            if (ce.subElements[FaceIndexes[f][(j+1)%3]]) {
                                newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, newEdgeCenterIt->second, face[(j+1)%3]});
                            }
                        } else if (ce.subElements[FaceIndexes[f][j]]) {
                            newTets.push_back(I4{tetCenterId, newFaceCenterIt->second, face[j], face[(j+1)%3]});
                        }
                    }
                } else if (ce.subElements[FaceIndexes[f][0]]) { // no face intersection, might have 0 or 1 edge cut
                    bool isSplit = false;
                    for (int j = 0; j < 3; ++j) {
                        auto sortedEdge = toI4<int,2>(sorted(std::array<int,2>{face[j], face[(j+1)%3]}));
                        const auto& newEdgeCenterIt = newNodeMapping.find(sortedEdge);
                        if (newEdgeCenterIt != newNodeMapping.end()) {
                            newTets.push_back(I4{tetCenterId, face[(j+2)%3], face[j], newEdgeCenterIt->second});
                            newTets.push_back(I4{tetCenterId, face[(j+2)%3], newEdgeCenterIt->second, face[(j+1)%3]});
                            isSplit = true;
                            break;
                        }
                    }
                    if (!isSplit) {
                        newTets.push_back(I4{tetCenterId, face[0], face[1], face[2]});
                    }
                }
            }
            newMesh[i] = newTets[0];
            for (int j = 1; j < newTets.size(); ++j) {
                newMesh.push_back(newTets[j]);
            }
        }
        return TetMesh<T>(move(newNodes), move(newMesh));
    }

public:
    TetMesh<T> tetMesh;
    TriMesh<T> triMesh;
    TetMesh<T> result;

    std::vector<std::array<T, 3>> *outNodes;
    std::vector<std::array<int, 4>> *outVoxels;
    std::vector<int> *outIDs;

    Cutter3D<T>(TetMesh<T> tetMesh, TriMesh<T> triMesh): tetMesh(tetMesh), triMesh(triMesh) {}
    Cutter3D<T>(std::vector<std::array<T,3>> volumeNodes, std::vector<std::array<int,4>> volumeVoxels,
                std::vector<std::array<T,3>> surfNodes,  std::vector<std::array<int,3>> surfFaces) {
        tetMesh.nodes_ = volumeNodes;
        tetMesh.mesh_ = volumeVoxels;
        tetMesh.initializeSurfaceMesh();
        tetMesh.computeConnectedComponents();

        triMesh.nodes_ = surfNodes;
        triMesh.mesh_ = surfFaces;
    }

    // static TetMesh<T> run(const TetMesh<T>& tetMesh, const TriMesh<T>& triMesh) {
    TetMesh<T> run() {
        TetBoundary2TetIds tetBoundary2TetIds;
        auto intersections = computeIntersections(tetMesh, triMesh, tetBoundary2TetIds);
        std::cout << "finished computing " << intersections.size() << " intersections\n";
//        for (auto& a: intersections) {
//            print<int,4>(a.first);
//            print<T,4>(a.second);
//        }
        std::set<int> cutTets;
        std::vector<CutElement> cutElements = split(tetMesh, intersections, tetBoundary2TetIds, cutTets);
//        for (auto& ce: cutElements) {
//            std::cout << ce.parentElementIndex << std::endl;
//            print<bool,4>(ce.subElements);
//        }
        std::vector<T3> newNodes;
        std::vector<I4> newMesh;
        merge(cutElements, tetMesh, newNodes, newMesh, intersections);
        std::cout << "finished split-merge\n";

        return subdivide(cutElements, tetMesh, newNodes, newMesh, intersections);
    }

    void runPlainOut() {
        TetBoundary2TetIds tetBoundary2TetIds;
        auto intersections = computeIntersections(tetMesh, triMesh, tetBoundary2TetIds);
        std::cout << "finished computing " << intersections.size() << " intersections\n";
//        for (auto& a: intersections) {
//            print<int,4>(a.first);
//            print<T,4>(a.second);
//        }
        std::set<int> cutTets;
        std::vector<CutElement> cutElements = split(tetMesh, intersections, tetBoundary2TetIds, cutTets);
//        for (auto& ce: cutElements) {
//            std::cout << ce.parentElementIndex << std::endl;
//            print<bool,4>(ce.subElements);
//        }
        std::vector<T3> newNodes;
        std::vector<I4> newMesh;
        merge(cutElements, tetMesh, newNodes, newMesh, intersections);
        std::cout << "finished split-merge\n";

        result = subdivide(cutElements, tetMesh, newNodes, newMesh, intersections);
        outNodes = &result.nodes_;
        outVoxels = &result.mesh_;
        outIDs = &result.connectedComponents_;
    }
};

#endif /* Cutting_h */
