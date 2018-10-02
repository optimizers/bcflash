%BND_LSQ_EXAMPLE Driver for running bcflash on random bounded least-squares
%
% Solve the problem
%    min_x |A*x - b|^2 subject to bL <= x <= bU
% for random data.
%
% PREREQUISITES:
%    Requires that directory containing +model/ is added to Matlab's path

clear all

% Set problem size
m = 2000; % number of rows
n = 500;  % number of variables
d = 0.1; % matrix density

% Define problem data
A  = sprand(m,n,d);
b  = randn(m,1);
bL = -rand(n,1);
bU =  rand(n,1);

% Build optimization problem
x0 = zeros(n,1);
nlp = model.nlpproblem('bnd_lsq', x0, zeros(0,1), zeros(0,1), bL, bU);

nlp.fobj_loc      = @(x) norm(A*x-b)^2;
nlp.gobj_loc      = @(x) A'*(A*x - b);
nlp.hlagprod_loc  = @(x,y,v) A'*(A*v);

% Create callback function (defined below)
callback = @(x,y,z,w) post_iteration(x,y,z,w);

% Set options
options.callback = callback;
options.maxiter  = 100;
options.gtolRel  = 1e-8;
options.gtolAbs  = 1e-8;
options.verbose  = 1;

% Build the solver
solver = bcflash(nlp, options);

% Run!
[x, info] = solver.solve(nlp.x0);

% Dummy post_iteration function
function [self, flag] = post_iteration(self, x, cgits, successful)
    flag = 0;
end