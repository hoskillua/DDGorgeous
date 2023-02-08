// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}
/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    if (!he.isInterior()) return 0;

    Vertex v0 = he.vertex();
    Vertex v1 = he.next().next().vertex();
    Vertex v2 = he.next().vertex();

    Vector3 x0 = inputVertexPositions[mesh.getVertexIndices()[v0]] - inputVertexPositions[mesh.getVertexIndices()[v1]];
    Vector3 x1 = inputVertexPositions[mesh.getVertexIndices()[v2]] - inputVertexPositions[mesh.getVertexIndices()[v1]];

    double D = dot(x1, x0);
    double C = cross(x1, x0).norm();

    if (C < 1E-9)
        C = 1E-9;

    return D/C;
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {


    /// Cotan Formula implementation was replaced for efficiency

    /*Halfedge iterh = v.halfedge();

    double S = 0;
    Vector3 xi = inputVertexPositions[mesh.getVertexIndices()[v]];


    do {

        Vector3 xj = inputVertexPositions[mesh.getVertexIndices()[iterh.next().vertex()]];
        Vector3 xk = inputVertexPositions[mesh.getVertexIndices()[iterh.next().next().vertex()]];

        S += (xk - xi).norm2() * cotan(iterh) + (xj - xi).norm2() * cotan(iterh.next().next());
        iterh = iterh.twin().next();
    } while (iterh != v.halfedge());

    return 0.125 * S;*/

    Halfedge iterh = v.halfedge();

    double S = 0;

    do {
        if(iterh.face().degree() == 3)
            S += faceArea(iterh.face());
        else
        {
            // area of polygons in xy plane
            Face f = iterh.face();
            Halfedge he = f.halfedge();
            double area = 0.0;

            // Calculate value of shoelace formula
            do
            {
                Vector3 pA = inputVertexPositions[he.vertex()];
                Vector3 pB = inputVertexPositions[he.next().vertex()];

                area += (pA.x + pB.x) * (pA.y + pB.y);

                he = he.next();
            } while (he != f.halfedge());

            // Return absolute value
            S += abs(area / 2.0);
        }
        iterh = iterh.twin().next();
    } while (iterh != v.halfedge() );

    return S/3.0;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    const Vector3 pi = inputVertexPositions[c.vertex()];
    const Vector3 pj = inputVertexPositions[c.halfedge().next().vertex()];
    const Vector3 pk = inputVertexPositions[c.halfedge().next().next().vertex()];

    const Vector3 vij = pj - pi;
    const Vector3 vik = pk - pi;

    return acos(dot(vij, vik) / (norm(vij) * norm(vik)));
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    // TODO
  const Vertex vi = he.tailVertex();
  const Vertex vj = he.tipVertex();

  Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[vi];

  const double eij_norm = norm(eij);

  eij /= eij_norm;

  const Face f1 = he.face();
  const Face f2 = he.twin().face();

  const Vector3 n1 = faceNormal(f1);
  const Vector3 n2 = faceNormal(f2);

  return atan2(dot(eij, cross(n1, n2)), dot(n1, n2));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
  // TODO


  Vector3 nAvg = { 0, 0, 0 };
  for (Face f : v.adjacentFaces())
    nAvg += faceNormal(f);

  nAvg = nAvg / norm(nAvg);

  return nAvg;
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    // TODO

  Vector3 nAvg = { 0, 0, 0 };
  for (Corner c : v.adjacentCorners())
    nAvg += angle(c) * faceNormal(c.face());

  nAvg = nAvg / norm(nAvg);

  return nAvg;
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO

  Vector3 nAvg = { 0, 0, 0 };
  for (Corner c : v.adjacentCorners())
  {
    const Vector3 pi = inputVertexPositions[c.vertex()];
    const Vector3 pj = inputVertexPositions[c.halfedge().next().vertex()];
    const Vector3 pk = inputVertexPositions[c.halfedge().next().next().vertex()];

    const Vector3 eij = pj - pi;
    const Vector3 eik = pk - pi;

    nAvg += cross(eij, eik) / (norm2(eij) * norm2(eik));
  }

  nAvg = nAvg / norm(nAvg);

  return nAvg;
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO

  Vector3 nAvg = { 0, 0, 0 };
  for (Face f : v.adjacentFaces())
    nAvg += faceArea(f) * faceNormal(f);

  nAvg = nAvg / norm(nAvg);

  return nAvg;
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
  Vector3 nGauss = { 0, 0, 0 };
  for (Halfedge he : v.incomingHalfedges())
  {
    // we are computing dihedral angle here instead of calling function to save edge vector calculations time
    const Vertex vi = he.tailVertex();
    const Vertex vj = he.tipVertex();

    Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[vi];

    const double eij_norm = norm(eij);

    eij /= eij_norm;

    const Face f1 = he.face();
    const Face f2 = he.twin().face();

    const Vector3 n1 = faceNormal(f1);
    const Vector3 n2 = faceNormal(f2);

    nGauss += atan2(dot(eij, cross(n1, n2)), dot(n1, n2)) * eij;
  }
  nGauss /= norm(nGauss);

  return nGauss;
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
  Vector3 nMean = { 0, 0, 0 };
  for (Halfedge he : v.incomingHalfedges())
  {
    // we are computing dihedral angle here instead of calling function to save edge vector calculations time
    const Vertex vi = he.tailVertex();
    const Vertex vj = he.tipVertex();

    Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[vi];


    nMean += (cotan(he) + cotan(he.twin())) * eij;
  }
  nMean /= norm(nMean);

  return nMean;
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
  double angleSum = 2 * PI;
  for (Corner c : v.adjacentCorners())
    angleSum -= angle(c);
  return angleSum;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
    // TODO
    double totalAngleDefect = 0;
    for (Vertex v : mesh.vertices())
    {
      totalAngleDefect += angleDefect(v);
    }
    return totalAngleDefect;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {

    // TODO
  double scalarMeanCurvature = 0;
  for (Halfedge he : v.incomingHalfedges())
  {
    // we are computing dihedral angle here instead of calling function to save edge vector calculations time
    const Vertex vi = he.tailVertex();
    const Vertex vj = he.tipVertex();

    Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[vi];

    const double eij_norm = norm(eij);

    scalarMeanCurvature += dihedralAngle(he) * eij_norm;
  }
  scalarMeanCurvature /= 2;
  return scalarMeanCurvature;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {

    // TODO

  double area = 0;
  for (Corner c : v.adjacentCorners())
  {
    const Vector3 pi = inputVertexPositions[c.vertex()];
    const Vector3 pj = inputVertexPositions[c.halfedge().next().vertex()];
    const Vector3 pk = inputVertexPositions[c.halfedge().next().next().vertex()];

    const Vector3 eij = pj - pi;
    const Vector3 eik = pk - pi;

    area += norm2(eij) * cotan(c.halfedge()) + norm2(eik) * cotan(c.halfedge().next().next());
  }
  area /= 8;

  return area;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
  double Hx2 = 2 * scalarMeanCurvature(v) / circumcentricDualArea(v);
  double G = angleDefect(v) / circumcentricDualArea(v);

  double k2 = (Hx2 + sqrt(Hx2 * Hx2 - 4 * G)) / 2;
  double k1 = Hx2 - k2;

  return std::make_pair(std::min(k1,k2), std::max(k1,k2)); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
  std::vector<Eigen::Triplet<double>> L_entries;

  for (Vertex v : mesh.vertices()) {
    double sum = 0;
    for (Halfedge he : v.incomingHalfedges())
    {
      double L_val = (cotan(he) + cotan(he.twin())) / 2;
      L_entries.push_back(Eigen::Triplet<double>(v.getIndex(), he.vertex().getIndex(), -L_val));
      sum += L_val;
    }
    L_entries.push_back(Eigen::Triplet<double>(v.getIndex(), v.getIndex(), sum + 1e-8));
  }

  //for (auto i : L_entries) std::cout << i.row() << " " << i.col() << " " << i.value() << "\n";

  Eigen::SparseMatrix<double> L(mesh.nVertices(), mesh.nVertices());
  L.setFromTriplets(L_entries.begin(), L_entries.end());
  return L;
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
  std::vector<Eigen::Triplet<double>> M_entries;

  for (Vertex v : mesh.vertices()) {
    M_entries.push_back(Eigen::Triplet<double>(v.getIndex(), v.getIndex(), barycentricDualArea(v)));
  }

  Eigen::SparseMatrix<double> M(mesh.nVertices(), mesh.nVertices());
  M.setFromTriplets(M_entries.begin(), M_entries.end());

  return M; // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
  std::vector<Eigen::Triplet<std::complex<double>>> Cl_entries;

  for (Vertex v : mesh.vertices()) {
    double sum = 0;
    for (Halfedge he : v.incomingHalfedges())
    {
      double L_val = (cotan(he) + cotan(he.twin())) / 2;
      std::complex<double> Cl_val(- L_val, 0);
      Cl_entries.push_back(Eigen::Triplet<std::complex<double>>(v.getIndex(), he.vertex().getIndex(), Cl_val));
      sum += L_val;
    }
    std::complex<double> Cl_sum(sum + 1e-8, 0);

    Cl_entries.push_back(Eigen::Triplet<std::complex<double>>(v.getIndex(), v.getIndex(), Cl_sum));
  }

  //for (auto i : L_entries) std::cout << i.row() << " " << i.col() << " " << i.value() << "\n";

  Eigen::SparseMatrix<std::complex<double>> Cl(mesh.nVertices(), mesh.nVertices());
  Cl.setFromTriplets(Cl_entries.begin(), Cl_entries.end());
  return Cl;

  return Cl; // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral