**luma**: library for unstructured adaptation
(c) Philip Claude Caplan, 2019

<img width="60px" src="doc/fig/luma.svg"/>

[![build status](https://gitlab.com/philipclaude/luma/badges/master/build.svg)](https://gitlab.com/philipclaude/luma/badges/master/build.svg)

**luma** is an unstructured mesh adaptation library with the following capabilities:

* dimension-independent mesh adaptation given a (1) mesh, (2) geometry description and (3) a metric field
* dimension-independent calculation of restricted Voronoi diagrams given (1) a set of sites and (2) a background mesh
* visualization of 2d, 3d and 4d meshes via (1) OpenGL and (2) websockets and WebGL

current/future developments include:
* parallel mesh adaptation using MPI
* optimal mesh to minimize the interpolation error in representing a function
* optimal mesh for representing a geometry (both static and time-dependent)
* dimension-independent calculation of quadrature rules using the centroidal Voronoi tessellation
* calculation of geometry-conforming Voronoi diagrams
* Python wrapper with an interface to blender
* curvilinear mesh adaptation

```
cd luma
mkdir build
mkdir build/release
cd build/release
cmake ../../
make luma
```

The libraries (**libluma.so** and **libluma.a** ) will be in **luma/build/release/lib**.

The main executable is **luma/build/release/bin/luma**.

Example 1: UGAWG Cube-Linear case
```
$ luma -adapt data/cube.mesh box Linear-3d tmp/cl.mesh
```

Example 2: UGAWG Cube-Cylinder Polar 2 case
```
$ luma -adapt ../data/cube-cylinder.mesh ../data/cube-cylinder.egads Polar2 ../tmp/ccp2.mesh
```

Example 3: mesh adaptation from a metric in a .solb file
```
luma -adapt input.mesh input.egads input.sol tmp/output.mesh
```

Example 4: visualization of a mesh (here, a 4d mesh from the Tesseract Wave case)
```
luma -plot wave.json
```

Notes:
* in ${EngSketchPad}/wvClient/WebViewer/wv-render.js, set 'preserveDrawingBuffer: true' so the canvas can be saved to png using toDataURL
* for linux: install xorg-dev if building the OpenGL visualizer
