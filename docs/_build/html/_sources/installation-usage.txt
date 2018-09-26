Installation and Basic Usage
============================

Dependencies
------------

BCFLASH uses the optimizers/model_ library to define optimization problems.

.. _model: https://github.com/optimizers/model

Installation
------------

To install BCFLASH, clone the the repository and add ``bcflash.m`` to your Matlab path.

Basic Usage
-----------

::

	nlp = model(...);       % model defining problem
	solver = bcflash(nlp);	% construct solver object.
	[x, info] = solver.solve(nlp.x0, varargin);	% begin solve

``nlpmodel`` contains many methods that are typically required by optimization solvers, such as objective function/gradient/Hessian evaluation, constraint/Jacobian/Hessian evaluation, Lagrangians, etc. However, BCFLASH requires that only the following functions to be implemented in ``model``:

* ``fobj(x)``: function value at current point, :math:`f(x)`.
* ``gobj(x)``: gradient at current point, :math:`\nabla f(x)`.
* ``hlagprod(x,y,v)``: Product with Lagrangian Hessian at current point, :math:`\nabla^2_x L(x,y)\cdot v`. Because the only constraints currently supported are bound constraints, this is equivalent to :math:`\nabla^2 f(x) \cdot v`.