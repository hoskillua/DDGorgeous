<div align="center">
<img src="./image/README/1676217625515.png" style="width:350px;" alt="Simplicial Complex">
</div>
<br />
<h1 align="center">DDGorgeous</h1>

This repo contains implementations of a collection of discrete differnetial geometry algorithms. They are based on a [C++ skeleton code](https://github.com/GeometryCollective/ddg-exercises) for the course assignments from [Discrete Differential Geometry](https://brickisland.net/DDGSpring2020/) (15-458/858).

This code framework uses [Geometry Central](https://github.com/nmwsharp/geometry-central) for geometry processing utilities and [Polyscope](https://github.com/nmwsharp/polyscope) for visualization, which were developed by Nick Sharp and others in the [Geometry Collective](http://geometry.cs.cmu.edu/). Also, It must be acknowledged that most of the illustrations used in this readme come from the course notes text provided with the mentioned course by Keenan Crane.

# Table of Contents

- [Results](#results)
  - [1.Simplicial Complex Operations](#1-simplicial-complex-operations)
  - [2.Discrete Exterior Calculus Operators](#2-discrete-exterior-calculus-operators)
  - [3.Normals & Curvatures](#3-normals--curvatures)
  - [4.The Laplace-Beltrami Operator & its Applications](#4-the-laplace-beltrami-operator--its-applications)
  - [5.Geodesics: The Heat Method](#5-geodesics-the-heat-method)
- [Dependencies](#dependencies-all-included)


# Results

Below are the highlights of the implemented algorithms.

## 1. Simplicial Complex Operations

<div align="center">
<img src="./image/README/Simplicial_complex_example.svg.png" style="width:350px;" alt="Simplicial Complex">
</div>
<br />

Given a mesh stored as a Halfedge structure, this part required building the incidence matrices and vector encodings for the mesh and its elements. These were then used to implement simple selection operations like Closure, Star, Link and boundary, as well as some simple boolen checks on a subcomplex of the mesh (`isComplex(()` and `isPureComplex()`)

<br />

|                                                        Operator                                                        |                                          Results (GIF)                                          |
| :---------------------------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------: |
|      **Star & Closure**<br />![img](image/README/1650891275357.png)<br />![img](image/README/1650891304489.png)      | Can use them together repeateadly to grow a selection.<br />![img](image/README/1650890503193.png) |
| **Link & Boundary**<br />![1676208491058](image/README/1676208491058.png)<br />![img](image/README/1650891245290.png) |  Can think of as an exclusive vs inclusive boundaries<br />![img](image/README/1650890508006.png)  |

## 2. Discrete Exterior Calculus Operators

<div align="center">
<img src="./image/README/1650891936774.png" style="width:450px;" alt="Discrete Exterior Calculus">
</div>
<br />

I got to learn about the discrete exterior calculus operators. The implementation of the discrete exterior derivative and the discrete Hodge star operators mainly involved implementing `cotan()` and `barycentricDualArea()` functions and using them to compute the discrete exterior derivative and the discrete Hodge star operators as matrices.

<br />

|                                                                     Operator                                                                     |            Results (GIF)            |
| :-----------------------------------------------------------------------------------------------------------------------------------------------: | :----------------------------------: |
| **Exterior Deravtive & Hodge Star**<br />![1676211127231](image/README/1676211127231.png)<br />![1676211109122](image/README/1676211109122.png) | ![img](image/README/1650902675076.png) |

## 3. Normals & Curvatures

<div align="center">
<img src="./image/README/1676211408743.png" alt="Normals & Curvatures">
</div>
<br />

I learned about a variety of ways to compute vertex normals, some of which are based on weighted averages of neighboring face normals. I also implemented the computation of the mean and Gaussian curvatures. These two were used to compute the principal curvatures.

<br />

|                                               Algorithm                                               |                                                                                 Results                                                                                 |
| :---------------------------------------------------------------------------------------------------: | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| **Vertex Normal Computation Methods<br /><br />![1676067070587](image/README/1676067070587.png)** |                                                              ![1676072347255](image/README/1676072347255.png)                                                              |
|     **Curvatures Computation <br /><br />![1676070772224](image/README/1676070772224.png)**     | **kmin & kmax:**<br />![1676070259394](image/README/1676070259394.png)<br />**Mean & Gaussian Curvature:**<br />![1676070266904](image/README/1676070266904.png) |

## 4. The Laplace-Beltrami Operator & its Applications

<div align="center">
<img src="./image/README/1676214620423.png" style="width:500px" alt="Laplace-Beltrami">
</div>
<br />

I got to implement the Laplace-Beltrami based on the cotangent Laplacian. It was used to implement the Poisson equation solver which was used to smoothely interpolate a function on the mesh. It was also used to implement mesh smoothing using the mean curvature flow and the stationary Laplacian mean curvature flow. The Complex Laplacian was used to implement the conformal parameterization of a mesh.

<br />

|                                              Algorithm                                              |                                                                                                                                                                                                                                              Results (GIF)                                                                                                                                                                                                                                              |
| :-------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|        **Poisson Equation<br /><br />![1676147668933](image/README/1676147668933.png)**        |                                                                                                                                                                                                                             ![1676148797658](image/README/1676148797658.png)                                                                                                                                                                                                                             |
| **Smoothing using Curvature Flows<br /><br />![1676147729135](image/README/1676147729135.png)** | **Mean Curvature Flow** (11 iterations)<br />(Updating Laplace Matrix in each iteration Vs using the initial one)<br />*Using the initial matrix (only updating mass matrix) helps with avoiding singularities*<br />![flow1](image/README/meancurvature1.gif)<br /><br />**Stationary-Laplacian Mean Curvature flow**<br />(~ 40 iterations, step size 0.001 vs 11 with step size 0.01)<br />*step size affects speed of convergance*<br />![flow2](image/README/meancurvature2.gif) |
|                                **Conformal Parameterization**                                |                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |

## 5. Geodesics: The Heat Method

<div align="center">
<img src="./image/README/1676214445987.png" style="width:350px" alt="The Heat Method">
</div>
<br />

Distance Computation was implemented using the heat method. This is based on a [paper by Crane et al](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/). First, the heat equation is solved on the mesh to compute the heat flow. Then, the heat flow is used to compute the distance from a given vertex to all other vertices. This involves normalizing the heat flow and negating it to get a vector field pointing along geodesics. A function whose gradient follows this vector field reproduces the final distance.

<br />

|                                                                                                                                                                                                                                                                     Algorithm                                                                                                                                                                                                                                                                     |                 Results (GIF)                 |
| :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :--------------------------------------------: |
| **Geodesics using**[ the Heat Method](https://www.cs.cmu.edu/~kmcrane/Projects/HeatMethod/) <br />![1676071513876](image/README/1676071513876.png)<br />![1676071566899](image/README/1676071566899.png)<br />(I) Heat is allowed to diffuse for short time (top-left). <br />(II) The temperature gradient (top-right) <br />is normalized & negated to get a field (bottom-left) pointing along geodesics. <br />(III) A function whose gradient follows the vector field recovers the final distance (bottom-right). | ![1676073261736](image/README/1676073261736.png) |

# Dependencies (all included)

1. Geometry processing and linear algebra - [Geometry Central](https://github.com/nmwsharp/geometry-central), which in turn has dependencies on [Eigen](https://eigen.tuxfamily.org) and/or [Suitesparse](https://people.engr.tamu.edu/davis/suitesparse.html).
2. Visualization - [Polyscope](https://github.com/nmwsharp/polyscope)
3. Unit tests - [Google Test](https://github.com/google/googletest)
