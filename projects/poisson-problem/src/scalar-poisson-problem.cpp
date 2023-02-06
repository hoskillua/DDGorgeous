// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"



/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    this->A = geometry->laplaceMatrix(); 
    this->M = geometry->massMatrix(); 
    this->totalArea = geometry->totalArea();
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {
    
    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
  double rhoBarValue = 0;
  for (int i = 0; i < rho.rows(); i++) {
    double area = geometry->barycentricDualArea((mesh->vertex(i)));
    rhoBarValue += rho[i] * area;
  }
  rhoBarValue = rhoBarValue / totalArea;
  Vector<double> rhoBar = rho;
  for (int i = 0; i < rho.rows(); i++) {
    rhoBar[i] = rhoBarValue;
  }
  Vector<double> rhs = -M * (rho - rhoBar);
  Eigen::SparseMatrix<double, 0, int> L = A;
  
  return solvePositiveDefinite(L, rhs);; // placeholder
}