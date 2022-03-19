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
        HV0[vi] = Eigen::Triplet<double>(vi, vi, barycentricDualArea(v));
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
        double length_ratio;
        if (e.isBoundary())
            if (e.halfedge().isInterior())
                length_ratio = cotan(e.halfedge().next().next());
            else
                length_ratio = cotan(e.halfedge().twin().next().next());
        else
            length_ratio = 0.5 * (cotan(e.halfedge().next().next()) + cotan(e.halfedge().twin().next().next()));
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

    std::vector<Eigen::Triplet<double>> EV1 (mesh.nEdges() * 2);
    size_t ei = 0, i = 0;
    for (Edge e : mesh.edges()) {
        ei = mesh.getEdgeIndices()[e];
        EV1[i] = Eigen::Triplet<double>( ei , mesh.getVertexIndices()[e.firstVertex()] , -1 );
        EV1[i + 1] = Eigen::Triplet<double>(ei, mesh.getVertexIndices()[e.secondVertex()], 1);
        i+=2;
    }
    Eigen::SparseMatrix<double> M0(mesh.nEdges(), mesh.nVertices());
    
    M0.setFromTriplets(EV1.begin(), EV1.end());

    return M0;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    std::vector<Eigen::Triplet<double>> FE1;
    size_t fi = 0;
    for (Face f : mesh.faces()) {
        fi = mesh.getFaceIndices()[f];
        Halfedge iterh = f.halfedge();
        for (int i = 0; i < f.degree(); i++)
        {
            FE1.push_back(Eigen::Triplet<double>(fi, mesh.getEdgeIndices()[iterh.edge()], -1 + (2 * (iterh.edge().firstVertex() == iterh.vertex()))));
            iterh = iterh.next();
        }
    }
    Eigen::SparseMatrix<double> M0(mesh.nFaces(), mesh.nEdges());

    M0.setFromTriplets(FE1.begin(), FE1.end());

    return M0;
}

} // namespace surface
} // namespace geometrycentral