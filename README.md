# CustomDists

Custom C++ samplers for various distributions that are missing from the C++ standard library and `boost/random`. 

## Generalize Inverse Gaussian (GIG) Distribution

The [GIGrvg](https://cran.r-project.org/web/packages/GIGrvg/index.html) R package implements a generalized inverse Gaussian sampler based on Hörmann and Leydold (2014), however:

* Using the `GIGrvg`'s low-level code in python or standalone C++ programs is challenging (one way to do this is to simply vendor code from the `GIGrvg` package into a new project with Python linkage, but this is still much more work than placing `LinkingTo: GIGrvg` in a `DESCRIPTION` file of an R package)
* The package / codebase is released under the GPL, making it incompatible with projects that use permissive licenses (MIT, BSD, Apache)

Scipy includes a [re-implementation of the Hörmann and Leydold (2014) algorithm](https://github.com/scipy/scipy/blob/92d2a8592782ee19a1161d0bf3fc2241ba78bb63/scipy/stats/_continuous_distns.py#L4933), released under the BSD license. The sampler is written primarily in pure python, with some Cython / C++ computation of special functions. 

We use this (nicely documented) python implementation as a blueprint for a C++ implementation which can be linked to R / Python or run in a standalone C++ program.
