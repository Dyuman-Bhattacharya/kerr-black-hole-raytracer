from setuptools import setup, Extension
import pybind11
import sys
import os

# Absolutely minimal include paths — no duplicates
include_dirs = [
    pybind11.get_include(),
    "cpp",
    os.path.join(sys.prefix, "include"),
]

ext = Extension(
    "kerr_cpp",
    sources=[
        "cpp/bindings.cpp",
        "cpp/kerr.cpp",
    ],
    include_dirs=include_dirs,
    language="c++",
    extra_compile_args=["/O2", "/openmp", "/std:c++17"],
    extra_link_args=["/openmp"],
)

setup(
    name="kerr_cpp",
    version="0.1",
    ext_modules=[ext],
)
