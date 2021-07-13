#!/bin/sh

git submodule init
git submodule update --remote src/third_party/pybind11
git submodule update src/third_party/webglpp
git submodule update --remote src/third_party/webglpp
