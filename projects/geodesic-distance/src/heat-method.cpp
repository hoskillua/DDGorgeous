// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "GeometryCentral/numerical/linear_solvers.h"
using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    double h = geometry->meanEdgeLength();
    h = h * h;

    // build the laplace matrix
    this->A = geometry->laplaceMatrix();
    // build the heat flow matrix
    this->F = geometry->massMatrix() + (h * A);
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {

    // TODO
    FaceData<Vector3>deltaU(*mesh, { 0, 0, 0 });
    
    for (Face f : mesh->faces()) {
      for (Halfedge he : f.adjacentHalfedges()) {
        deltaU[f] += u[he.next().tipVertex().getIndex()] * cross(geometry->faceNormal(f), geometry->halfedgeVector(he));
      }
      deltaU[f] = - deltaU[f].normalize();
    }
    return deltaU; // placeholder
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {

    // TODO
    Vector<double> divX = Vector<double>::Zero(mesh->nVertices()); // placeholder
    // Note: core/geometry.cpp has faceArea() function
    // Note: core/geometry.cpp has cotan() function
    // compute the divergence of X as ∇·X = 0.5 * sum over j cotθ1(e1 · Xj) +cotθ2(e2 · Xj)
    // we can loop over all edges and compute the contribution of each edge to the divergence of X
    for(Edge e : mesh->edges())
    {
        // get the two faces adjacent to the edge
        Halfedge he = e.halfedge();
        Halfedge he_twin = he.twin();
        Face f1 = he.face();
        Face f2 = he_twin.face();
        // get the cotangent of the angles at the two vertices of the edge
        double cot_theta1 = geometry->cotan(he);
        double cot_theta2 = geometry->cotan(he_twin);
        // get the indices of the two vertices of the edge
        int v1 = he.vertex().getIndex();
        int v2 = he_twin.vertex().getIndex();
        // get the vector X at the two faces
        Vector3 X1 = X[f1];
        Vector3 X2 = X[f2];
        // get the edge vector
        Vector3 e_vec = geometry->halfedgeVector(he);
        // compute the contribution of the edge to the divergence of X
        double divV1 = 0.5 * (cot_theta1 * dot(e_vec, X1) + cot_theta2 * dot(e_vec, X2));
        double divV2 = 0.5 * (cot_theta1 * dot(- e_vec, X1) + cot_theta2 * dot(- e_vec, X2));
        // add the contribution to the divergence of X
        divX[v1] += divV1;
        divX[v2] += divV2;
    }
    return divX;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // TODO
    // Solve the heat equation ∂u/∂t = Δu
    // use u(x, t) = u0(x) + ∫_0^t e^(-Δt) * u(x, 0) dt
    // since u(x, 0) = 0, we have u(x, t) = ∫_0^t e^(-Δt) * u0(x) dt
    // we can solve this equation by solving the linear system F * u = delta
    // where F = M + h * A, M is the mass matrix, A is the laplace matrix, and h is the mean edge length
    Vector<double> ut = solvePositiveDefinite(SparseMatrix<double>(F), delta);
    // compute the vector field X = -∇u / |∇u|
    FaceData<Vector3> X = computeVectorField(ut);
    // compute the divergence of X
    Vector<double> divX = computeDivergence(X);
    // compute the geodesic distance Δφ = ∇.X
    Vector<double> phi = - solvePositiveDefinite(SparseMatrix<double>(A), divX);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}