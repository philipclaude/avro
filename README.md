**avro**
(c) Philip Claude Caplan, 2017-2021

<img width="60px" src="doc/fig/avro.svg"/>

[![build status](https://gitlab.com/philipclaude/avro/badges/main/pipeline.svg)](https://gitlab.com/philipclaude/avro/badges/main/pipeline.svg)

[![coverage](https://gitlab.com/philip/avro/badges/main/coverage.svg)](https://gitlab.com/philip/avro/badges/main/coverage.svg)

**avro** is an unstructured mesh adaptation library with the following capabilities:

* dimension-independent parallel mesh adaptation given a (1) mesh, (2) geometry description and (3) a metric field
* dimension-independent calculation of restricted Voronoi diagrams given (1) a set of sites and (2) a background mesh
* visualization of 2d, 3d and 4d meshes via (1) OpenGL and (2) websockets and WebGL

documentation: https://philipclaude.gitlab.io/avro/

quickstart:

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

Example 3: visualization of a mesh
```
avro -plot ../data/cube-cylinder.mesh
```

Notes:
* for linux: install xorg-dev if building the OpenGL visualizer
