**ursa**: unstructured adaptation library
(c) Philip Claude Caplan, 2019

<img width="60px" src="doc/fig/ursa.png"/>

[![build status](https://gitlab.com/philipclaude/ursa/badges/master/build.svg)](https://gitlab.com/philipclaude/ursa/badges/master/build.svg)

**ursa** is an unstructured mesh adaptation library with the following capabilities:

* dimension-independent mesh adaptation given a (1) mesh, (2) geometry description and (3) a metric field
* dimension-independent calculation of restricted Voronoi diagrams given (1) a set of sites and (2) a background mesh
* visualization of 2d, 3d and 4d meshes in the browser via (1) OpenGL and (2) websockets and WebGL

```
cd ursa
mkdir build
mkdir build/release
cd build/release
cmake ../../
make ursa
```

The libraries (**libursa.so** and **libursa.a** ) will be in **ursa/build/release/lib**.

The main executable is **ursa/build/release/bin/ursa**.

Example 1: UGAWG Cube-Linear case
```
$ ursa -adapt data/cube.mesh box Linear-3d tmp/cl.mesh
```

Example 2: UGAWG Cube-Cylinder Polar 2 case
```
$ ursa -adapt ../data/cube-cylinder.mesh ../data/cube-cylinder.egads Polar2 ../tmp/ccp2.mesh
```

Example 3: mesh adaptation from a metric in a .solb file
```
ursa -adapt input.mesh input.egads input.sol tmp/output.mesh
```

Example 4: visualization of a mesh (here, a 4d mesh from the Tesseract Wave case)
```
ursa -plot wave.json
```

Notes:
* in ${EngSketchPad}/wvClient/WebViewer/wv-render.js, set 'preserveDrawingBuffer: true' so the canvas can be saved to png using toDataURL
* for linux: install xorg-dev if building the OpenGL visualizer
