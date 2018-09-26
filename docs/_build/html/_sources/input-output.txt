Input/Output
============
Arguments are passed into the solver constructor. Optional argumentas are passed in as Matlab's ``varargin`` using name-value pairs. For example:
::

	solver = bcflash(nlp, 'maxiter', 100, 'gtolRel', 1e-8);
	[x, info] = solver.solve(x0);


Input
-----

============== ======== ================ ===========
Name           Required Default          Description
============== ======== ================ ===========
x0             required                  initial point
maxiter	       optional :math:`10n`      max number of outer iterations
maxcgiter      optional :math:`n`        max CG iterations for subproblem
cgtol          optional :math:`0.1`      subproblem residual tolerance
fatol          optional :math:`0`        min allowable abs. function reduction
frtol          optional :math:`10^{-12}` min allowable rel. function reduction
min_radius     optional :math:`10^{-16}` min trust-region radius
gtolRel        optional :math:`10^{-6}`  relative gradient tolerance
gtolAbs        optional :math:`10^{-6}`  absolute gradient tolerance
stop_tol       optional []               if empty, overwritten by gtolRel and gtolAbs
fmin           optional :math:`-10^{32}` min allowable function value
mu0            optional :math:`0.01`     sufficient decrease parameter
verbose        optional 1                display log?
fid            optional 1                output file
callback       optional []               see below
exit_user_only optional 0                see below
============== ======== ================ ===========

* callback: A callback function that gets called at the end of every iteration, which allows custom monitoring and logging, as well as the ability to modify the underlying problem mid-solve. It has the function signature ::

	[self, flag] = callback(self, x, cgits, successful)

  where the inputs are:

  * x: the current iterate,
  * cgits: the total number of CG iterations so far,
  * successful: boolean flag indicating if trust-region step was accepted or rejected.

  The output (other than bcflash itself) is ``flag``, where if:

  * ``flag=0``: solver continues uninterrupted,
  * ``flag=7``: solver immediately terminates,
  * ``flag=8``: the underlying optimization has been modified. The solver takes this into account and continues.

* exit_user_only: Ignore eFlag=1,2, so that termination occurs only if callback returns 1, or an error occurs.

Output
------

BCFLASH returns the latest iterate, ``x``, and a struct, ``info``, with the following fields:

======= ===========
Field   Description
======= ===========
eFlag   exit flag (see below)
msg     output message
obj     final objective value
gpnorm  norm of projected gradient
iters   total number of outer iterations
cgiters total number of CG iterations
======= ===========

Exit Codes
----------

BCFLASH terminates when one of the following conditions is met:

* eFlag = 1 (Optimal solution found)
	:math:`\|P(\nabla f(x_k)\| \le \mbox{gtolAbs} + \|\nabla f(x_0)\| \cdot \mbox{gtolRel}`
* eFlag = 3 (Unbounded below)
	:math:`f(x_k) \le \mbox{fmin}`
* eFlag = 4 (Absolute function tolerance)
	predicted or actual reduction less than fatol
* eFlag = 5 (Relative function tolerance)
	predicted or actual reduction less than :math:`\mbox{frtol} \cdot |f(x_k)|`
* eFlag = 6 (Trust region radius too small)
	trust region radius less than min_radius
* eFlag = 2 (Too many iterations)
	exceeded number of outer iterations
* eFlag = 7 (User requested exit)
	callback returned ``flag=1``