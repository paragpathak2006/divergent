// divergent_test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include<vector>
#include<unordered_set>
#include<unordered_map>
#include "Input_output/Loader.h"
// ---------------------------------------------------------------------------------
// Implement a simple mesh data structure capable of representing an all - triangle mesh and its constituent vertices.
// This data structure should include functions to do the following :

//Geometry
struct Point {
    double x, y, z;
};
inline double dist(Point P, Point Q) {
    return (P.x - Q.x) * (P.x - Q.x) + (P.y - Q.y) * (P.y - Q.y) + (P.z - Q.z) * (P.z - Q.z);
}


//Topology
class Vertex {
public:
    // d. Return the coordinates of a given vertex.

    Point point;
    // b. Given a vertex , return the adjacent faces / vertices.
    std::vector<int> adjVertices;
    std::vector<int> adjFaces;

    bool is_deleted = false;
    // e. Delete a vertex, with optional flag to delete all connected faces(if a vertex).
    Vertex(double X, double Y, double Z) {
        point.x = X;
        point.y = Y;
        point.z = Z;
        is_deleted = false;
    }

};
inline double dist(Vertex P, Vertex Q) {return dist(P.point, Q.point);}

typedef std::vector<Vertex> Vertices;

class Face {
public:

    // b. Given a face, return the adjacent faces / vertices.
    int vertices[3];
    std::vector<int> adjFaces;
    bool is_deleted = false;


    // e. Delete a face
    void remove() {
        is_deleted = true;
    }


    // f. Construct a new face from vertices, and a new vertex from coordinates
    Face(int v0, int v1, int v2) {
        vertices[0] = v0;
        vertices[1] = v1;
        vertices[2] = v2;
        is_deleted = false;
    }


    // g. Flip the sense of a face.
    void flip_sense() {
        std::swap(vertices[1], vertices[2]);
    }


    // 4. Write a function that returns all faces with minimum angle below a specified angle in degrees.
    bool check_below_threshold_angle(const Vertices& v, const double& cos2_angle) const {
        const Vertex& i = v[vertices[0]];
        const Vertex& j = v[vertices[1]];
        const Vertex& k = v[vertices[2]];

        const double a = dist(i, j);
        const double b = dist(j, k);
        const double c = dist(k, i);

        if (c + b > a) { double cos2_a = (c + b - a) * (c + b - a) / (4 * c * b); if (cos2_angle > cos2_a) return true; }
        if (a + c > b) { double cos2_b = (a + c - b) * (a + c - b) / (4 * a * c); if (cos2_angle > cos2_b) return true; }
        if (a + b > c) { double cos2_c = (a + b - c) * (a + b - c) / (4 * a * b); if (cos2_angle > cos2_c) return true; }

        return false;

    }
};
typedef std::vector<Face> Faces;


void remove(Vertex v, Faces& faces, bool flag = true) {
    v.is_deleted = true;
    if (flag)
        for (auto f : v.adjFaces)
            faces[f].remove();
}



class Mesh {
// c. Return all the vertices or faces.
public:
    Vertices vertices;
    Faces faces;
};

class Edge {
public:
    int i, j;
    Edge(int ii, int jj) { i = ii; j = jj; }
    bool operator==(const Edge& b) const { return i == b.i && j == b.j; }
};

// Edge Hash functor
struct Hash_of_Edge {
    size_t operator() (const Edge& P) const {
        return P.i * 18397 + P.j * 20483;
    }
};

typedef std::unordered_set<Edge, Hash_of_Edge> Edges;

// ---------------------------------------------------------------------------------
// 2. Write a function that returns whether all faces are consistently oriented.
bool checkOriention(Edges& edges, int i, int j) {
    Edge edge(i, j);

    if (edges.count(edge) != 0)
        return false;

    edges.emplace(i,j);
    return true;
}

bool checkOriention(const Faces& faces) {
    std::unordered_set<Edge, Hash_of_Edge> edges;
    edges.reserve(faces.size() * 3);
    for (const auto& face : faces) {
        if (!checkOriention(edges, face.vertices[0], face.vertices[1])) return false;
        if (!checkOriention(edges, face.vertices[1], face.vertices[2])) return false;
        if (!checkOriention(edges, face.vertices[2], face.vertices[0])) return false;
    }
}

// ---------------------------------------------------------------------------------
// 3. Write a function that returns the number of loops bounding a surface mesh
void addEdges(Edges& edges, int i, int j,int& num_edges) {
    Edge edge(j, i);

    if (edges.count(edge) != 0)
        return;

    edges.emplace(i, j);
    num_edges++;
}

typedef std::unordered_multimap<int, Edge> VertexToEdge;

int num_loops(const Faces& faces, int num_vertices, int num_faces) {
    std::unordered_set<Edge, Hash_of_Edge> edges;
    int num_edges = 0;
    edges.reserve(faces.size() * 3);
    for (const auto& face : faces) {
        addEdges(edges, face.vertices[0], face.vertices[1], num_edges);
        addEdges(edges, face.vertices[1], face.vertices[2], num_edges);
        addEdges(edges, face.vertices[2], face.vertices[0], num_edges);
    }

    int euler_characteristic = num_vertices - num_edges + num_faces;


    return 2 - euler_characteristic;
    //VertexToEdge vertexToEdgemap;
    //vertexToEdgemap.reserve(edges.size() * 2);
    //for (const auto& edge : edges) {
    //    vertexToEdgemap.emplace(edge.i, edge);
    //    vertexToEdgemap.emplace(edge.j, edge);
    //}

    //auto i = edges.begin()->i;
    //auto j = edges.begin()->j;
}


// ---------------------------------------------------------------------------------
// 5. Write a function that collapses all edges with length below a specified threshold.
void Edge_Colapse(Faces faces, Vertices v, int i,int j){
    for (auto f : v[j].adjFaces) 
        for(int ii = 0; ii<3;ii++) 
            if (faces[f].vertices[ii] == j) 
                faces[f].vertices[ii] = i;
    
    for (const auto& adjV : v[j].adjVertices)
        for (auto& adjV2 : v[adjV].adjVertices)
            if (adjV2 == j) 
                adjV2 = i;

    // add adjFaces to i
}


// ---------------------------------------------------------------------------------
// 6. (stretch)Write a function that performs a diagonal edge swap for triangles having an obtuse angle and an angle below a specified threshold in degrees.
void get_opposite_vertices(int i, int& j, int& k) {
    if (i == 0) { j == 1; k == 2; return; }
    if (i == 1) { j == 2; k == 0; return; }
    if (i == 2) { j == 0; k == 1; return; }
}

void diagonal_edge_swap(Face& f1, Face& f2) {
    bool found = false;
    int j1, k1;
    int j2, k2;
    for (size_t i1 = 0; i1 < 3; i1++)
    {
        get_opposite_vertices(i1, j1, k1);
        for (size_t i2 = 0; i2 < 3; i2++)
        {
            get_opposite_vertices(i2, j2, k2);
            if (f1.vertices[j1] == f2.vertices[k2] && f1.vertices[k1] == f2.vertices[j2])
            {
                f1.vertices[k1] = f2.vertices[i2];
                f2.vertices[k2] = f1.vertices[i1];
                return;
            }
        }

    }

}


// ---------------------------------------------------------------------------------
// a. Read/write basic .obj files; sample files are included with this specification. Code can assume mesh is all triangles. Code should be capable of reading files with smoothing group ('s'), object names ('o'), polygon groups ('g'), material libraries ('mtllib'), and use materials ('usemtl'), but does not need to provide query access to those data, nor write those data.

void get_points(const objl::Mesh& curMesh, Vertices& points) {
    int length = curMesh.Vertices.size();

    points.reserve(length);
    for (int j = 0; j < length; j++) {
        points.emplace_back(
            curMesh.Vertices[j].Position.X,
            curMesh.Vertices[j].Position.Y,
            curMesh.Vertices[j].Position.Z
        );
    }
}

void get_faces(const objl::Mesh& curMesh, Faces& faces) {
    int length = curMesh.Vertices.size();
    faces.reserve(length);
    for (int j = 0; j < length; j += 3)
        faces.emplace_back(curMesh.Indices[j], curMesh.Indices[j + 1], curMesh.Indices[j + 2]);
}

// ---------------------------------------------------------------------------------

int main()
{
    objl::Mesh mesh;        get_mesh("shuttle.obj", mesh);
    Vertices vertices;      get_points(mesh, vertices);
    Faces faces;            get_faces(mesh, faces);
    
    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
