**avro**

(c) Philip Claude Caplan, 2017-2021

<img width="60px" src="doc/fig/avro.svg"/>

[![build status](https://gitlab.com/philipclaude/avro/badges/main/pipeline.svg)](https://gitlab.com/philipclaude/avro/badges/main/pipeline.svg)

**avro** is an unstructured mesh adaptation library with the following capabilities:

* parallel mesh adaptation in 2d, 3d and 4d, given a (1) mesh, (2) geometry description and (3) a metric field,
* dimension-independent calculation of restricted Voronoi diagrams given (1) a set of sites and (2) a background mesh,
* dimension-independent computation of a semi-discrete optimal transport map using a Newton-based method,
* visualization of high-order solution fields on curved 2d, 3d and linear 4d meshes.

documentation: https://philipclaude.gitlab.io/avro/
coverage: https://philipclaude.gitlab.io/avro/coverage_results/

quickstart (after installing dependencies described in documentation):

```
git clone --recursive https://gitlab.com/philipclaude/avro.git
cd avro
mkdir build
mkdir build/release
cd build/release
cmake ../../
make avro
```

To compile avro without any external dependencies (some functionalities will be disabled), set the cmake variable avro_LITE = true.

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

For any bugs or feature requests, please submit an issue in this GitLab project.

If you use **avro**, I'd really like to hear from you!
Please send me an email at pcaplan@middlebury.edu.
