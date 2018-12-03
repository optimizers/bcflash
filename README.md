# Bound-Constrained Flash (BCFLASH)

## Overview

BCFLASH is an optimization solver for

`` min_x f(x) subject to l <= x <= u``

that is essentially a Matlab implementation of <a href="https://www.mcs.anl.gov/~more/tron/">TRON</a>.

## Documentation

Hosted on <a href="https://bcflash.readthedocs.io/en/latest/index.html">ReadTheDocs</a>.

## Installation

BCFLASH uses the <a href="https://github.com/optimizers/model">optimizers/model</a> to define the optimization problems.

Once the model library is installed, then installing BCFLASH simply involves cloning this repository and adding bcflash.m to your path.

## Basic Usage

```Matlab
nlp = model(...); % model defining problem
solver = bcflash(nlp); % construct solver object.
[x, info] = solver.solve(nlp.x0); % begin solve
```
