# CamparyBindings

This started as a binding of
[Campary](https://homepages.laas.fr/mmjoldes/campary/), which is released under GPLv2.

With time it became a port of Campary. It is not production ready but many of the algorithms from the Ph. D. thesis of Valentina Popescu
[Towards fast and certified multiple-precision libraries](https://hal.archives-ouvertes.fr/tel-01534090v2) are implemented.

An important disclaimer is that I am not related in any form to the original authors so I take the blame for anything not working as it should.

## The rationale behind the library
The library implements multiprecision through the use of vectors
of Floating point numbers.
This is implemented in classical libraries as the double-double library or the qd library.

The idea behind Campary is to use vectors of double, error free transformations and GPU to implement a fast multiprecision library.

## Basic usage
It is possible to construct a CamparyFloat by

```julia
x = CamparyFloat{4}(1.0)
y = CamparyFloat{4}([1.0, 2^-54, 2^-108, 0])
```

there is a lot of syntactic sugar, so most operations work out of the box.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://orkolorko.github.io/CamparyBindings.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://orkolorko.github.io/CamparyBindings.jl/dev)
[![Build Status](https://travis-ci.com/orkolorko/CamparyBindings.jl.svg?branch=master)](https://travis-ci.com/orkolorko/CamparyBindings.jl)
