# DDGorgeous

This repo contains implementations of a collection of discrete differnetial geometry algorithms. They are based on a [C++ skeleton code](https://github.com/GeometryCollective/ddg-exercises) for course assignments from [Discrete Differential Geometry](https://brickisland.net/DDGSpring2020/) (15-458/858).

This code framework uses [Geometry Central](https://github.com/nmwsharp/geometry-central) for geometry processing utilities and [Polyscope](https://github.com/nmwsharp/polyscope) for visualization, which were developed by Nick Sharp and others in the [Geometry Collective](http://geometry.cs.cmu.edu/). Also, It must be acknowledged that most of the illustrations used in this readme are from the notes paper provided with the mentioned course by Keenan Crane.

## Results

| Algorithms                                                                                 |                                                                       Illustration                                                                       | Result (GIF)                                                                                                    |
| :----------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------: | --------------------------------------------------------------------------------------------------------------- |
| Simplicial<br />Complex <br />Operations <br />(Star, <br />Closure, <br />Link, Boundary) | ![img](image/README/1650888559097.png) <br />![img](image/README/1650888578915.png) ![img](image/README/1650888604129.png) ![img](image/README/1650888865111.png) | ![](image/README/1650890398226.png)<br />![](image/README/1650890503193.png)<br />![](image/README/1650890508006.png) |

## Dependencies (all included)

1. Geometry processing and linear algebra - [Geometry Central](https://github.com/nmwsharp/geometry-central), which in turn has dependencies on [Eigen](https://eigen.tuxfamily.org) and/or [Suitesparse](https://people.engr.tamu.edu/davis/suitesparse.html).
2. Visualization - [Polyscope](https://github.com/nmwsharp/polyscope)
3. Unit tests - [Google Test](https://github.com/google/googletest)
