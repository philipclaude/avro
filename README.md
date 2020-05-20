**avro**: library for unstructured mesh adaptation
(c) Philip Claude Caplan, 2019-2020

<img width="60px" src="doc/fig/avro.svg"/>


[![build status](https://gitlab.com/philipclaude/avro/badges/master/pipeline.svg)](https://gitlab.com/philipclaude/avro/badges/master/pipeline.svg)

[![coverage](https://gitlab.com/philipclaude/avro/badges/master/coverage.svg)](https://gitlab.com/philipclaude/avro/badges/master/coverage.svg)

**avro** is an unstructured mesh adaptation library with the following capabilities:

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
cd avro
mkdir build
mkdir build/release
cd build/release
cmake ../../
make avro
```

The libraries (**libavro.so** and **libavro.a** ) will be in **avro/build/release/lib**.

The main executable is **avro/build/release/bin/avro**.

Example 1: UGAWG Cube-Linear case
```
$ avro -adapt data/cube.mesh box Linear-3d tmp/cl.mesh
```

Example 2: UGAWG Cube-Cylinder Polar 2 case
```
$ avro -adapt ../data/cube-cylinder.mesh ../data/cube-cylinder.egads Polar2 ../tmp/ccp2.mesh
```

Example 3: mesh adaptation from a metric in a .solb file
```
avro -adapt input.mesh input.egads input.sol tmp/output.mesh
```

Example 4: visualization of a mesh (here, a 4d mesh from the Tesseract Wave case)
```
avro -plot wave.json
```

Notes:
* in ${EngSketchPad}/wvClient/WebViewer/wv-render.js, set 'preserveDrawingBuffer: true' so the canvas can be saved to png using toDataURL
* for linux: install xorg-dev if building the OpenGL visualizer
