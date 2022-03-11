// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    /// TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    std::vector<Eigen::Triplet<size_t>> EV1 (mesh->nEdges() * 2);
    size_t ei = 0, i = 0;
    for (Edge e : mesh->edges()) {
        ei = mesh->getEdgeIndices()[e];
        EV1[i] = Eigen::Triplet<size_t>( ei , mesh->getVertexIndices()[e.firstVertex()] , 1 );
        EV1[i + 1] = Eigen::Triplet<size_t>(ei, mesh->getVertexIndices()[e.secondVertex()], 1);
        i+=2;
    }
    Eigen::SparseMatrix<size_t> M0(mesh->nEdges(), mesh->nVertices());
    
    M0.setFromTriplets(EV1.begin(), EV1.end());

    return M0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    std::vector<Eigen::Triplet<size_t>> FE1;
    size_t fi = 0;
    for (Face f : mesh->faces()) {
        fi = mesh->getFaceIndices()[f];
        Halfedge iterh = f.halfedge();
        for (int i = 0; i < f.degree(); i++)
        {
            FE1.push_back(Eigen::Triplet<size_t>(fi, mesh->getEdgeIndices()[iterh.edge()], 1));
            iterh = iterh.next();
        }
    }
    Eigen::SparseMatrix<size_t> M0(mesh->nFaces(), mesh->nEdges());

    M0.setFromTriplets(FE1.begin(), FE1.end());

    return M0; 
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    Vector <size_t> V;
    V.resize(mesh->nVertices(), 1);
    for (Vertex v : mesh->vertices())
        V[mesh->getVertexIndices()[v]] = 0;
    for (size_t v : subset.vertices)
        V[v] = 1;
    return V;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    Vector <size_t> E;
    E.resize(mesh->nEdges(), 1);
    for (Edge e : mesh->edges())
        E[mesh->getEdgeIndices()[e]] = 0;
    for (size_t e : subset.edges)
        E[e] = 1;
    return E;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {


    Vector <size_t> F;
    F.resize(mesh->nFaces(), 1);
    for (Face f : mesh->faces())
        F[mesh->getFaceIndices()[f]] = 0;
    for (size_t f : subset.faces)
        F[f] = 1;
    return F;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    MeshSubset newsubset = subset.deepCopy();

    Vector<size_t> E0 = A0 * buildVertexVector(subset);
    for (Edge e : mesh->edges()) 
    {
        int i = mesh->getEdgeIndices()[e];
        if (E0[i]) newsubset.edges.insert(i);
    }
    Vector<size_t> F0 = A1 * buildEdgeVector(newsubset);
    for (Face f : mesh->faces())
    {
        int i = mesh->getFaceIndices()[f];
        if (F0[i]) newsubset.faces.insert(i);
    }

    return newsubset; 
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    MeshSubset newsubset = subset.deepCopy();

    Vector<size_t> E0 = A1.transpose() * buildFaceVector(subset);
    for (Edge e : mesh->edges())
    {
        int i = mesh->getEdgeIndices()[e];
        if (E0[i]) newsubset.edges.insert(i);
    }
    Vector<size_t> V0 = A0.transpose() * buildEdgeVector(newsubset);
    for (Vertex v : mesh->vertices())
    {
        int i = mesh->getVertexIndices()[v];
        if (V0[i]) newsubset.vertices.insert(i);
    }

    return newsubset; 
    
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset out = closure(star(subset));
    MeshSubset in = star(closure(subset));
    
    for (size_t v : in.vertices)
        out.deleteVertex(v);
    for (size_t e : in.edges)
        out.deleteEdge(e);
    for (size_t f : in.faces)
        out.deleteFace(f);

    return out;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    /// TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    /// TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    MeshSubset newsubset;
    int deg = //isPureComplex(subset);
        1;
    if (deg == 2)
    {
        Vector<size_t> E0 = A1.transpose() * buildFaceVector(subset);
        for (Edge e : mesh->edges())
        {
            int i = mesh->getEdgeIndices()[e];
            if (E0[i] == 1) newsubset.edges.insert(i);
        }
        return closure(newsubset);
    }
    else if (deg == 1)
    {
        Vector<size_t> V0 = A0.transpose() * buildEdgeVector(subset);
        for (Vertex v : mesh->vertices())
        {
            int i = mesh->getVertexIndices()[v];
            if (V0[i] == 1) newsubset.vertices.insert(i);
        }
        return closure(newsubset);
    }
    else if (deg == 0)
        return newsubset;
    else
    {
        /// error message ?
        return subset;
    }
}