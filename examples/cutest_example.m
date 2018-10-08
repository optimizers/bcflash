%CUTEST_EXAMPLE Driver for running bcflash on CUTEst test problems
%
% PREREQUISITES:
%    Requires that cutest-mirror be installed from here:
%       https://github.com/optimizers/cutest-mirror
%    and that CUTEst problems are compiled
%
%    Requires that directory containing +model/ is added to Matlab's path

clear all

% Path to directory containing compiled CUTEst problems
% For example, if problem 3pk is to be run, then the compiled problem
% should be in directory: cutest_path/3pk
cutest_path = '';
% Problem name
pname = 'ssc';

fname = [cutest_path '/' pname];

% Build optimization problem
nlp = model.cutestmodel(fname, true);

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