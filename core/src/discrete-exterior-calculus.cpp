// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    std::vector<Eigen::Triplet<double>> HV0(mesh.nVertices());

    for (Vertex v : mesh.vertices())
    {
        int vi = mesh.getVertexIndices()[v];
        HV0[vi] = Eigen::Triplet<double>(vi, vi, v.isBoundary()? 1 : barycentricDualArea(v));;
    }

    Eigen::SparseMatrix<double> H0(mesh.nVertices(), mesh.nVertices());

    H0.setFromTriplets(HV0.begin(), HV0.end()); 

    return H0;
    //return identityMatrix<double>(1); // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    std::vector<Eigen::Triplet<double>> HV1(mesh.nEdges());

    for (Edge e : mesh.edges())
    {
        int ei = mesh.getEdgeIndices()[e];
        double length_ratio = e.isBoundary() ? 1 : 0.5 * (cotan(e.halfedge().next().next()) + cotan(e.halfedge().twin().next().next()));
        HV1[ei] = Eigen::Triplet<double>(ei, ei, length_ratio);
    }

    Eigen::SparseMatrix<double> H1(mesh.nEdges(), mesh.nEdges());

    H1.setFromTriplets(HV1.begin(), HV1.end());

    return H1;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    std::vector<Eigen::Triplet<double>> HV2(mesh.nFaces());

    for (Face f : mesh.faces())
    {
        int fi = mesh.getFaceIndices()[f];
        double area = faceArea(f);
        if (area < 1e-9)
            area = 1e-9;
        //double area = 0.5 * (cotan(f.halfedge().next().next()) + cotan(f.halfedge().twin().next().next()));
        HV2[fi] = Eigen::Triplet<double>(fi, fi, 1/area);
    }

    Eigen::SparseMatrix<double> H2(mesh.nFaces(), mesh.nFaces());

    H2.setFromTriplets(HV2.begin(), HV2.end());

    return H2;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

} // namespace surface
} // namespace geometrycentral