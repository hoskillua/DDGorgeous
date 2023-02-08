// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {
    
    // TODO
  SparseMatrix<std::complex<double>> Ec = 0.5 * geometry->complexLaplaceMatrix();

  std::vector<Eigen::Triplet<std::complex<double>>> A_entries;
  
  for(BoundaryLoop bl : mesh->boundaryLoops())
    for (Halfedge he : bl.adjacentHalfedges())
    {
      A_entries.push_back(
        Eigen::Triplet<std::complex<double>>(
          he.tailVertex().getIndex(),
          he.tipVertex().getIndex(), 
          std::complex<double>(0, 0.25)
          )
      );
      A_entries.push_back(
        Eigen::Triplet<std::complex<double>>(
          he.tipVertex().getIndex(),
          he.tailVertex().getIndex(),
          std::complex<double>(0, -0.25)
          )
      );
    }

  Eigen::SparseMatrix<std::complex<double>> A(Ec.rows(), Ec.cols());
  A.setFromTriplets(A_entries.begin(), A_entries.end());

  Ec -= A;

  return Ec; // placeholder
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    // TODO
    SparseMatrix<std::complex<double>> Ec = buildConformalEnergy();

    // Solve the eigenvalue problem with the smallest eigenvalue using solveInversePowerMethod
    Vector<std::complex<double>> eig_vec;
    eig_vec = solveInversePowerMethod(Ec);

    VertexData<Vector2> result(*mesh);

    for (Vertex v: mesh->vertices())
    {
      std::complex<double> result_cmplx = eig_vec[v.getIndex()];
      result[v].x = result_cmplx.real();
      result[v].y = result_cmplx.imag();
    }

    return result; // placeholder
}