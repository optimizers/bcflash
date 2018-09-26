classdef bcflash
   
% Copyright (C) 2018  Michael P. Friedlander, Dominique Orban, and Ron Estrin
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

   properties (Access = protected)
      x              % solution
      fx             % objective value at x
      pgnorm         % norm of projected gradient at x
      time_total     % total solve time
      eFlag          % exit flag
      exit_msg       % string indicating exit
      nlp            % copy of the nlp object
      mu0            % sufficient decrease parameter
      stop_tol       % stopping tolerance
      gtolRel        % relative tolerance in projected gradient
      gtolAbs        % absolute tolerance in projected gradient
      cgtol          % CG tolerance for subproblems
      fatol          % absoulte error in function
      frtol          % relative error in function
      fmin           % min allowed function value for declaring unbounded
      fid            % File ID of where to direct log output
      verbose        % log level
      gnorm0         % norm of the gradient at x0
      num_successful_itns = 0 % number of successful iterations
      iteration = 0  % iteration counter
      cgiters = 0    % total number of CG iterations
      maxiter        % maximum number of iterations
      maxcgiter      % maximum number of CG iterations per Newton step
      min_radius     % Minimum trust region radius
      callback       % function called at the end of each iteration
      exit_user_only % If true, solver stops only if post_iteration returns 1
   end
   
   properties (Access = private, Constant)
      EXIT_NONE                  = 0;
      EXIT_OPTIMAL               = 1;
      EXIT_ITERATIONS            = 2;
      EXIT_UNBOUNDED             = 3;
      EXIT_FATOL                 = 4;
      EXIT_FRTOL                 = 5;
      EXIT_MIN_RADIUS            = 6;
      EXIT_USER_REQUEST          = 7;
      EXIT_RESTART               = 8;
      EXIT_UNKNOWN               = 9;
      EXIT_MSG = {
         'Optimal solution found'
         'Too many iterations'
         'Unbounded below'
         'Absolute function tolerance'
         'Relative function tolerance'
         'Trust region radius too small'
         'User requested exit'
         'User requested restart'
         'Unknown exit'};
      
      % Constants used to manipulate the TR radius. These are the numbers
      % used by TRON.
      sig1 = 0.25;
      sig2 = 0.50;
      sig3 = 4.00;
      eta0 = 1e-4;
      eta1 = 0.25;
      eta2 = 0.75;
      
      % Log header and body formats.
      logH = '\n%5s%1s %13s  %13s  %5s(%3s)  %10s  %10s\n';
      logB = '%5i%1s %13.6e  %13.6e  %5i(%3s)  %10.3e  %10.3e  %3s\n';
      logT = {'iter','','f(x)','|g(x)|','cg','ext','preRed','radius'};
      
   end
   
   methods (Access = public)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = bcflash(nlp, varargin)
         % ---------------------------------------------------------------------
         % Parse input parameters and initialize local variables.
         % ---------------------------------------------------------------------
         p = inputParser;
         p.addParameter('maxiter', 10*length(nlp.x0));
         p.addParameter('maxcgiter', length(nlp.x0));
         p.addParameter('cgtol', 0.1);
         p.addParameter('fatol', 0);
         p.addParameter('frtol', 1e-12);
         p.addParameter('min_radius', 1e-16);
         p.addParameter('gtolRel', 1e-6);
         p.addParameter('gtolAbs', 1e-6);
         p.addParameter('stop_tol', []);
         p.addParameter('fmin', -1e32);
         p.addParameter('mu0', 0.01);
         p.addParameter('verbose', 1);
         p.addParameter('fid', 1);
         p.addParameter('exit_user_only', false);
         p.addParameter('callback', ...
             @(x,y,z,w) bcflash.post_iteration(x,y,z,w), ...
             @(f) isa(f, 'function_handle'));
         p.parse(varargin{:});

         % ---------------------------------------------------------------------
         % Store various objects and parameters.
         % ---------------------------------------------------------------------
         self.nlp = nlp;
         self.gtolRel = p.Results.gtolRel;
         self.gtolAbs = p.Results.gtolAbs;
         self.cgtol = p.Results.cgtol;
         self.fatol = p.Results.fatol;
         self.frtol = p.Results.frtol;
         self.min_radius = p.Results.min_radius;
         self.stop_tol = p.Results.stop_tol;
         self.maxiter = p.Results.maxiter;
         self.maxcgiter = p.Results.maxcgiter;
         self.fmin = p.Results.fmin;
         self.mu0 = p.Results.mu0;
         self.verbose = p.Results.verbose;
         self.fid = p.Results.fid;
         self.exit_user_only = p.Results.exit_user_only;
         self.callback = p.Results.callback;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, info] = solve(self, x)

         [self, x, f, g, delta] = self.pre_solve(x);

         xl = self.nlp.bL;
         xu = self.nlp.bU;
         
         % Actual and predicted reductions. Initial inf value prevents
         % exits based on related on first iteration.
         actred = inf;
         prered = inf;
         
         % Miscellaneous iteration
         alphac = 1;
         sigma1 = self.sig1;
         sigma2 = self.sig2;
         sigma3 = self.sig3;
         self.eFlag = self.EXIT_NONE;
         successful = [];
         cgits = 0;
         cginfo = '';
         
         while true
            
            % ---------------------------------------------------------
            % Check stopping conditions.
            % ---------------------------------------------------------

            % Optimality. This test takes precedence, and so overrides
            % other previously set exit conditions.
            exit = ~self.exit_user_only && (self.pgnorm <= self.stop_tol);
            if exit
               self.eFlag = self.EXIT_OPTIMAL;
            end

            % Unboundedness.
            exit = f < self.fmin;
            if ~self.eFlag && exit
               self.eFlag = self.EXIT_UNBOUNDED;
            end

            % Actual and predicted absolute reductions small.
            exit = max( abs(actred), prered) <= self.fatol;
            if ~self.eFlag && exit
               self.eFlag = self.EXIT_FATOL;
            end
            
            % Actual and predicted relative reductions small. 
            exit = max( abs(actred), prered) <= self.frtol*abs(f);
            if ~self.eFlag && exit
               self.eFlag = self.EXIT_FRTOL;
            end
            
            % Trust region radius is too small
            exit = delta <= self.min_radius;
            if ~self.eFlag && exit
               self.eFlag = self.EXIT_MIN_RADIUS; 
            end
            
            % Iterations.
            exit = ~self.exit_user_only && (self.iteration >= self.maxiter);
            if ~self.eFlag && exit
               self.eFlag = self.EXIT_ITERATIONS;
            end
            
            % Restart flag. This resets so as not to exit.
            if self.eFlag == self.EXIT_RESTART
               self.eFlag = self.EXIT_NONE;
               flag = '*';
            else
               flag = '';
            end
            
            % ---------------------------------------------------------
            % Print current iteration to log.
            % ---------------------------------------------------------
            if self.verbose
               if mod(self.iteration, 20) == 0
                  self.printf(self.logH, self.logT{:});
               end
               accrej = '';
               if ~successful
                  accrej = 'rej';
               end
               self.printf(self.logB, self.iteration, flag, f, ...
                  self.pgnorm, cgits, cginfo, prered, delta, accrej);
            end

            % ---------------------------------------------------------
            % Act on exit conditions.
            % ---------------------------------------------------------
            if self.eFlag
               self.exit_msg = self.EXIT_MSG{self.eFlag};
               self.x = x;
               self.fx = f;
               break
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Iteration starts here.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            self.iteration = self.iteration + 1;
            
            xc = x;
            fc = f;
            gc = g;
            
            % Hessian operator at current iterate.
            Aprod = @(v)self.nlp.hlagprod(x, zeros(0,1), v);
            
            % Cauchy step: s.
            [alphac, s] = self.cauchy(x, gc, Aprod, delta, alphac);
            
            % Projected Newton iteration: x.
            [x, s, cgits, cginfo] = self.spcg(Aprod, x, gc, delta, ...
               self.cgtol, s, self.maxcgiter, xl, xu);
            self.cgiters = self.cgiters + cgits;
            
            % Predicted reduction.
            As = self.nlp.hlagprod(x, zeros(0,1), s);
            prered = -(s'*gc + 0.5*s'*As);
            
            % Compute the objective at the new x.
            f = self.nlp.fobj(x);
            actred = fc - f;
            snorm = norm(s);
            
            % Update the trust-region radius.
            if self.num_successful_itns == 0
               delta = min(delta, snorm);
            end
            gts = gc'*s;
            if f - fc - gts <= 0
               alpha = sigma3;
            else
               alpha = max(sigma1, -0.5*gts/(f-fc-gts));
            end
            
            successful = actred >= self.eta0*prered;

            if actred < self.eta0*prered
               delta = min(max(alpha, sigma1)*snorm, sigma2*delta);
            elseif actred < self.eta1*prered
               delta = max(sigma1*delta, min(alpha*snorm, sigma2*delta));
            elseif actred < self.eta2*prered
               delta = max(sigma1*delta, min(alpha*snorm, sigma3*delta));
            else
               delta = max(delta, min(alpha*snorm, sigma3*delta));
            end
            
            if successful
               % Successful iteration. Evaluate the gradient at the new x.
               self.num_successful_itns = self.num_successful_itns + 1;
               g = self.nlp.gobj(x);
               self.pgnorm = bcflash.gpnrm2(x, xl, xu, g);
            end               
                        
            % ---------------------------------------------------------
            % Post-iteration function.
            % ---------------------------------------------------------
            [self, cflag] = self.callback(self, x, cgits, successful);
            if cflag == self.EXIT_RESTART
               % The user may have redefined the objective function.
               % Re-evaluate the objective and gradient.
               if successful
                  f = self.nlp.fobj(x);
                  g = self.nlp.gobj(x);
                  self.pgnorm = bcflash.gpnrm2(x, xl, xu, g);
               else
                  fc = self.nlp.fobj(xc);
                  gc = self.nlp.gobj(xc);
               end
               self.eFlag = cflag;
            elseif cflag == self.EXIT_USER_REQUEST
               self.eFlag =  self.EXIT_USER_REQUEST;
            end
            
            if ~successful
               x = xc;
               f = fc;
               g = gc;
            end            
            
         end % loop
         
         self = self.post_solve();
         
         info.eFlag  = self.eFlag;
         info.msg    = self.exit_msg;
         info.obj    = self.fx;
         info.gpnorm = self.pgnorm;
         info.iters  = self.iteration;
         info.cgiter = self.cgiters;
         
      end % function solve
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function val = getprop(self, prop)
         %GETVAL  Get a property value from this object.
         % Needed because the properties for this class are private.
         val = self.(prop);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = setprop(self, prop, val)
         %SETVAL  Get a property value from this object.
         % Needed because the properties for this class are private.
         self.(prop) = val;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end % public methods
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   methods (Access = private)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [self, x, f, g, delta] = pre_solve(self, x)
         
         xl = self.nlp.bL;
         xu = self.nlp.bU;
         
         % Make sure initial point is feasible before evaluating function.
         x = bcflash.mid(x, xl, xu);
         
         % First objective and gradient evaluation.
         f = self.nlp.fobj(x);
         g = self.nlp.gobj(x);
         
         % Initial trust-region radius.
         gnorm = norm(g);
         delta = gnorm;

         % Initialize stopping tolerance.
         self.pgnorm = bcflash.gpnrm2(x, xl, xu, g);
         self.gnorm0 = gnorm;
         if ~isempty(self.stop_tol)
            % User overrode the default. Don't reset it.
         else
            self.stop_tol = self.gtolAbs + self.gtolRel * gnorm;
         end
         
         % ---------------------------------------------------------------------
         % Print header.
         % ---------------------------------------------------------------------
         if self.verbose
            self.printf('\n');
            self.printf('%s\n',repmat('=',1,80));
            self.printf('Bound-Constrained FLASH \n');
            self.printf('%s\n\n',repmat('=',1,80));
            self.printf(self.nlp.formatting())
            self.printf('\nParameters\n----------\n')
            self.printf('%-15s: %3s %8i'  ,'iter max','',self.maxiter);
            self.printf('%5s','');
            self.printf('%-15s: %3s %8.1e\n','mu0','',self.mu0);
            self.printf('%-15s: %3s %8.1e'  ,'frtol','',self.frtol);
            self.printf('%5s','');
            self.printf('%-15s: %3s %8.1e\n'  ,'fmin','',self.fmin);
            self.printf('%-15s: %3s %8.1e'  ,'cgtol','',self.cgtol);
            self.printf('%5s','');
            self.printf('%-15s: %3s %8.1e\n'  ,'rel stop tol','',self.gtolRel);
            self.printf('%-15s: %3s %8.1e'    ,'fatol','',self.fatol);
            self.printf('%5s','');
            self.printf('%-15s: %3s %8.1e\n','abs stop tol','',self.gtolAbs);
            self.printf('%-15s: %3s %8s'  ,'','','');
            self.printf('%5s','');
            self.printf('%-15s: %3s %8.1e\n'  ,'stopping tol','',self.stop_tol);
            self.printf('\n');
         end
         
         % ---------------------------------------------------------------------
         % Start the clock.
         % ---------------------------------------------------------------------
         self.time_total = tic;

      end % function pre_solve
      
      function self = post_solve(self)
         
         % ---------------------------------------------------------------------
         % Stop the clock.
         % ---------------------------------------------------------------------
         self.time_total = toc(self.time_total);
         
         % ---------------------------------------------------------------------
         % Print footer.
         % ---------------------------------------------------------------------
         if self.verbose
            self.printf('\n EXIT: %s\n', self.exit_msg);
            self.printf('\n')
            self.printf(' %-27s  %6i     %-17s  %15.8e\n',...
               'No. of iterations', self.iteration,...
               'Objective value', self.fx);
            t1 = self.nlp.ncalls_fobj + self.nlp.ncalls_fcon;
            t2 = self.nlp.ncalls_gobj + self.nlp.ncalls_gcon;
            self.printf(' %-27s  %6i     %-17s    %6i\n',...
               'No. of calls to objective' , t1,...
               'No. of calls to gradient', t2);
            self.printf(' %-27s  %6i     %-22s  %10.2e\n',...
               'No. of Hessian-vector prods', self.nlp.ncalls_hvp,...
               'No. of successful iterations', self.num_successful_itns);
            self.printf('\n');
            tt = self.time_total;
            t1 = self.nlp.time_fobj+self.nlp.time_fcon; t1t = round(100*t1/tt);
            t2 = self.nlp.time_gcon+self.nlp.time_gcon; t2t = round(100*t2/tt);
            self.printf(' %-24s %6.2f (%3d%%)  %-20s %6.2f (%3d%%)\n',...
               'Time: function evals' , t1, t1t,...
               'gradient evals', t2, t2t);
            t1 = self.nlp.time_hvp; t1t = round(100*t1/tt);
            self.printf(' %-24s %6.2f (%3d%%)  %-20s %6.2f (%3d%%)\n',...
               'Time: Hessian-vec prods', t1, t1t,...
               'total solve', tt, 100);
         end
                  
      end % function post_solve
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [alpha, s] = cauchy(self, x, g, Aprod, delta, alpha)
         %CAUCHY
         %
         % Computes a Cauchy step that satisfies a trust-region constraint
         % and a sufficient-decrease condition.
         %
         % The Cauchy step is computed for the quadratic
         %
         %       q(s) = 0.5*s'*A*s + g'*s,
         %
         % where A is a symmetric operator, and g is a vector. Given a
         % parameter alpha, the Cauchy step is
         %
         %       s[alpha] = P[x - alpha*g] - x,
         %
         % with P the projection onto the n-dimensional interval [xl,xu].
         % The Cauchy step satisfies the trust-region constraint and the
         % sufficient-decrease condition
         %
         %       |s| <= delta,      q(s) <= mu_0*(g'*s),
         %
         % where mu_0 is a constant in (0,1).

         interpf =  0.1;     % interpolation factor
         extrapf = 10.0;     % extrapolation factor
         xl = self.nlp.bL;
         xu = self.nlp.bU;
         
         % Find the minimal and maximal break-point on x - alpha*g.
         [~, ~, brptmax] = bcflash.breakpt(x, -g, xl, xu);
         
         % Evaluate the initial alpha and decide if the algorithm
         % must interpolate or extrapolate.
         s = bcflash.gpstep(x, -alpha, g, xl, xu);
         
         if norm(s) > delta
            interp = true;
         else
            As = Aprod(s);
            gts = g'*s;
            q = 0.5*s'*As + gts;
            interp = (q >= self.mu0*gts);
         end
         
         % Either interpolate or extrapolate to find a successful step.
         if interp
            
            % Reduce alpha until a successful step is found.
            search = true;
            while search
               
               alpha = interpf*alpha;               
               s = bcflash.gpstep(x, -alpha, g, xl, xu);
               if norm(s) <= delta
                  As = Aprod(s);
                  gts = g'*s;
                  q = 0.5*s'*As + gts;
                  search = (q >= self.mu0*gts);
               end
            end
            
         else
            
            % Increase alpha until a successful step is found.
            search = true;
            alphas = alpha;
            while search && alpha <= brptmax
               
               alpha = extrapf*alpha;
               s = bcflash.gpstep(x, -alpha, g, xl, xu);
               if norm(s) <= delta
                  As = Aprod(s);
                  gts = g'*s;
                  q = 0.5*s'*As + gts;
                  if q <= self.mu0*gts
                     search = true;
                     alphas = alpha;
                  end
               else
                  search = false;
               end
            end
            
            % Recover the last successful step.
            alpha = alphas;
            s = bcflash.gpstep(x, -alpha, g, xl, xu);
         end
         
      end % function cauchy
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, w] = prsrch(self, Aprod, x, g, w, xl, xu)
         %PRSRCH  Projected search.
         %
         % [x, w] = prsrch(Aprod, x, g, w, xl, xu) where
         %
         %     Inputs:
         %     Aprod is a function handle to compute matrix-vector products
         %     x        current point
         %     g        current gradient
         %     w        search direction
         %     xl       vector of lower bounds
         %     xu       vector of upper bounds
         %     mu0      linesearch parameter
         %     interpf  interpolation parameter
         %
         %     Output:
         %     x is the final point P[x + alpha*w]
         %     w is the step s[alpha]
         %
         %     This subroutine uses a projected search to compute a step
         %     that satisfies a sufficient decrease condition for the quadratic
         %
         %           q(s) = 0.5*s'*A*s + g'*s,
         %
         %     where A is a symmetric matrix in compressed column storage,
         %     and g is a vector. Given the parameter alpha, the step is
         %
         %           s[alpha] = P[x + alpha*w] - x,
         %
         %     where w is the search direction and P the projection onto the
         %     n-dimensional interval [xl,xu]. The final step s = s[alpha]
         %     satisfies the sufficient decrease condition
         %
         %           q(s) <= mu_0*(g'*s),
         %
         %     where mu_0 is a constant in (0,1).
         %
         %     The search direction w must be a descent direction for the
         %     quadratic q at x such that the quadratic is decreasing
         %     in the ray  x + alpha*w for 0 <= alpha <= 1.

         interpf = 0.5; % Interpolation factor
         
         % Set the initial alpha = 1 because the quadratic function is
         % decreasing in the ray x + alpha*w for 0 <= alpha <= 1.
         alpha = 1;
         nsteps = 0;
         
         % Find the smallest break-point on the ray x + alpha*w.
         [~, brptmin, ~] = bcflash.breakpt(x, w, xl, xu);
         
         % Reduce alpha until the sufficient decrease condition is
         % satisfied or x + alpha*w is feasible.
         search = true;
         while search && alpha > brptmin
            
            % Calculate P[x + alpha*w] - x and check the sufficient
            % decrease condition.
            nsteps = nsteps + 1;
            s = bcflash.gpstep(x, alpha, w, xl, xu);
            As = Aprod(s);
            gts = g'*s;
            q = 0.5*s'*As + gts;
            
            if q <= self.mu0*gts
               search = false;
               
            else
               % This is a crude interpolation procedure that
               % will be replaced in future versions of the code.
               alpha = interpf*alpha;
               
            end
         end
         
         % Force at least one more constraint to be added to the active
         % set if alpha < brptmin and the full step is not successful.
         % There is sufficient decrease because the quadratic function
         % is decreasing in the ray x + alpha*w for 0 <= alpha <= 1.
         if alpha < 1 && alpha < brptmin
            alpha = brptmin;
         end
         
         % Compute the final iterate and step.
         s = bcflash.gpstep(x, alpha, w, xl, xu);
         x = bcflash.mid(x + alpha*w, xl, xu);
         w = s;
         
      end % function prsrch
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [x, s, iters, info] = spcg(self, Aprod, x, g, delta, rtol, s, itermax, xl, xu)
         %SPCG  Minimize a bound-constraint quadratic.
         %
         % This subroutine generates a sequence of approximate minimizers
         % for the subproblem
         %
         %       min { q(x) : xl <= x <= xu }.
         %
         % The quadratic is defined by
         %
         %       q(x[0]+s) = 0.5*s'*A*s + g'*s,
         %
         % where x[0] is a base point provided by the user, A is a
         % symmetric operator, and g is a vector.
         %
         % At each stage we have an approximate minimizer x[k], and generate
         % a direction p[k] by using a preconditioned conjugate gradient
         % method on the subproblem
         %
         %       min { q(x[k]+p) : | L'*p | <= delta, s(fixed) = 0 },
         %
         % where fixed is the set of variables fixed at x[k], delta is the
         % trust-region bound, and L is an incomplete Cholesky factorization
         % of the submatrix
         %
         %       B = A(free:free),
         %
         % where free is the set of free variables at x[k]. Given p[k],
         % the next minimizer x[k+1] is generated by a projected search.
         %
         % The starting point for this subroutine is x[1] = x[0] + s, where
         % x[0] is a base point and s is the Cauchy step.
         %
         % The subroutine converges when the step s satisfies
         %
         %       | (g + A*s)[free] | <= rtol*| g[free] |
         %
         % In this case the final x is an approximate minimizer in the
         % face defined by the free variables.
         %
         % The subroutine terminates when the trust region-bound does not
         % allow further progress, that is, |L'*p[k]| = delta. In this
         % case, the final x satisfies q(x) < q(x[k]).
         %
         % On exit info is set as follows:
         %
         %      info = 'cnv'  Convergence. The final step s satisfies
         %                |(g + A*s)[free]| <= rtol*|g[free]|,
         %                and the final x is an approximate minimizer
         %                in the face defined by the free variables.
         %
         %      info = 'bnd'  Termination. The trust-region bound does
         %                not allow further progress.
         %
         %      info = 'max'  Failure to converge within itermax iterations.
         
         n = length(x);
         
         % Compute A*(x[1] - x[0]) and store in w.
         As = Aprod(s);
         
         % Store current iterate as center of trust-region center
         x0 = x;
         
         % Compute the Cauchy point.
         x = bcflash.mid(x + s, xl, xu);
         
         % Start the main iteration loop.
         % There are at most n iterations because at each iteration
         % at least one variable becomes active.
         iters = 0;
         info = 'max';
         for nfaces = 1:n
            
            % Determine the free variables at the current minimizer.
            % The indices of the free variables are stored in the first
            % n free positions of the array indfree.
            indfree = (xl < x) & (x < xu);
            nfree = sum(indfree);
            
            % Exit if there are no free constraints.
            if nfree == 0
               info = 'cnv';
               return
            end
            
            % Compute the gradient grad q(x[k]) = g + A*(x[k] - x[0]),
            % of q at x[k] for the free variables.
            % Recall that w contains  A*(x[k] - x[0]).
            % Compute the norm of the reduced gradient Z'*g.
            wa = g(indfree);
            gfree = As(indfree) + wa;
            gfnorm = norm(wa);
            
            % Solve the trust-region subproblem in the free variables
            % to generate a direction p[k]. Store p[k] in the array w.
            tol = rtol*gfnorm;
            stol = 0;
            
            % Create the submatrix operator.
            Bprod = @(x)bcflash.Afree(x, Aprod, indfree, n);
            
            % Fix delta to keep iterates within trust region for trpcg
            d = x-x0;
            dd = d(~indfree);
            delta_new = sqrt(delta^2 - dd'*dd);
            dd = d(indfree);
            
            L = speye(nfree); % No preconditioner for now.
            [w, itertr, infotr] = bcflash.trpcg(dd, Bprod, gfree, delta_new, ...
               L, tol, stol, itermax);
            
            iters = iters + itertr;
            w = L' \ w;
            
            % Use a projected search to obtain the next iterate.
            % The projected search algorithm stores s[k] in w.
            xfree = x(indfree);
            xlfree = xl(indfree);
            xufree = xu(indfree);
            
            [xfree, w] = self.prsrch(Bprod, xfree, gfree, w, xlfree, xufree);
            
            % Update the minimizer and the step.
            % Note that s now contains x[k+1] - x[0].
            x(indfree) = xfree;
            s(indfree) = s(indfree) + w;
            
            % Compute A*(x[k+1] - x[0]) and store in w.
            As = Aprod(s);
            
            % Compute the gradient grad q(x[k+1]) = g + A*(x[k+1] - x[0])
            % of q at x[k+1] for the free variables.
            gfree = As(indfree) + g(indfree);
            gfnormf = norm(gfree);
            
            % Convergence and termination test.
            % We terminate if the preconditioned conjugate gradient method
            % encounters a direction of negative curvature, or
            % if the step is at the trust region bound.
            
            if gfnormf <= rtol*gfnorm
               info = 'cnv';
               return
            elseif infotr == 3 || infotr == 4
               info = 'bnd';
               return
            elseif iters > itermax
               info = 'max';
               return
            end
            
         end % for faces
         
      end % function spcg
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function printf(self, varargin)
         fprintf(self.fid, varargin{:});
      end % function printf

   end % private methods
  
   methods(Static)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function xFree = Afree(xFree, Aprod, indfree, n)
         z = zeros(n,1);
         z(indfree) = xFree;
         z = Aprod(z);
         xFree = z(indfree);
      end % function Afree
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [w, iters, info] = trpcg(d, Aprod, g, delta, L, tol, stol, itermax)
         %TRPCG
         %
         % This subroutine uses a preconditioned conjugate gradient method
         % to find an approximate minimizer of the trust region subproblem
         %
         %       min { q(s) : || L'*s || <= delta }.
         %
         % where q is the quadratic
         %
         %       q(s) = 0.5*s'*A*s + g'*s,
         %
         % A is a symmetric operator, L is a lower triangular matrix, and g
         % is a vector.
         %
         % This subroutine generates the conjugate gradient iterates for
         % the equivalent problem
         %
         %       min { Q(w) : || w || <= delta }.
         %
         % where Q is the quadratic defined by
         %
         %       Q(w) = q(s),      w = L'*s.
         %
         % Termination occurs if the conjugate gradient iterates leave
         % the trust region, a negative curvature direction is generated,
         % or one of the following two convergence tests is satisfied.
         %
         % Convergence in the original variables:
         %
         %       || grad q(s) || <= tol
         %
         % Convergence in the scaled variables:
         %
         %       || grad Q(w) || <= stol
         %
         % Note that if w = L'*s, then L*grad Q(w) = grad q(s).
         %
         % On exit info is set as follows:
         %
         %       info = 1  Convergence in the original variables.
         %                 || grad q(s) || <= tol
         %
         %       info = 2  Convergence in the scaled variables.
         %                 || grad Q(w) || <= stol
         %
         %       info = 3  Negative curvature direction generated.
         %                 In this case || w || = delta and a direction
         %                 of negative curvature w can be recovered by
         %                 solving L'*w = p.
         %
         %       info = 4  Conjugate gradient iterates exit the
         %                 trust region. In this case || w || = delta.
         %
         %       info = 5  Failure to converge within itermax iterations.

         n = length(g);
         
         % Initialize the iterate w and the residual r.
         w = zeros(n,1);
         
         % Initialize the residual t of grad q to -g.
         % Initialize the residual r of grad Q by solving L*r = -g.
         % Note that t = L*r.
         t = -g;
         r = L \ -g;
         
         % Initialize the direction p.
         p = r;
         
         % Initialize rho and the norms of r and t.
         rho = r'*r;
         rnorm0 = sqrt(rho);
         
         % Exit if g = 0.
         if rnorm0 == 0
            iters = 0;
            info = 1;
            return
         end
         
         for iters = 1:itermax
            
            % Compute z by solving L'*z = p.
            z = L' \ p;
            
            % Compute q by solving L*q = A*z and save L*q for
            % use in updating the residual t.
            Az = Aprod(z);
            z = Az;
            q = L \ Az;
            
            % Compute alpha and determine sigma such that the trust region
            % constraint || w + sigma*p || = delta is satisfied.
            ptq = p'*q;
            if ptq > 0
               alpha = rho/ptq;
            else
               alpha = 0;
            end
            sigma = bcflash.trqsol(d+w, p, delta);
            
            % Exit if there is negative curvature or if the
            % iterates exit the trust region.
            if (ptq <= 0 || alpha >= sigma)
               if sigma ~= 0
                  w = w + sigma*p;
               end
               if ptq <= 0
                  info = 3;
               else
                  info = 4;
               end
               
               return
               
            end
            
            % Update w and the residuals r and t.
            % Note that t = L*r.
            w = w + alpha*p;
            r = r - alpha*q;
            t = t - alpha*z;
            
            % Exit if the residual convergence test is satisfied.
            rtr = r'*r;
            rnorm = sqrt(rtr);
            tnorm = norm(t);
            
            if tnorm <= tol
               info = 1;
               return
            end
            
            if rnorm <= stol
               info = 2;
               return
            end
            
            % Compute p = r + beta*p and update rho.
            beta = rtr/rho;
            p = r + beta*p;
            rho = rtr;
            
         end
         
         info = 5;
         
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function x = mid(x, xl, xu)
         %MID  Project a vector onto the box defined by [xl, xu].
         x = max( x, xl );
         x = min( x, xu );
      end % function mid
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function s = gpstep(x, alpha, w, xl, xu)
         %GPSTEP  Compute the gradient projection step.
         %
         % Compute the gradient-projection step
         %
         % s = P[x + alpha*w] - x,
         %
         % where P is the projection onto the box [xl, xu].

         aw = alpha*w;
         s = x + aw;
         
         iLow = s < xl;         % violate lower bound
         iUpp = s > xu;         % violate upper bound
         iFre = ~(iLow | iUpp); % free
         
         s(iLow) = xl(iLow) - x(iLow);
         s(iUpp) = xu(iUpp) - x(iUpp);
         s(iFre) = aw(iFre);
                  
      end % function gpstep
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function pnorm = gpnrm2(x, xl, xu, g)

         nfix = xl < xu;
         low = nfix & x == xl;
         upp = nfix & x == xu;
         fre = nfix & ~(low | upp);
         
         pnorm1 = norm(g(low & g < 0))^2;
         pnorm2 = norm(g(upp & g > 0))^2;
         pnorm3 = norm(g(fre))^2;
         
         pnorm = sqrt(pnorm1 + pnorm2 + pnorm3);
         
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function [nbrpt, brptmin, brptmax] = breakpt(x, w, xl, xu)
         %BREAKPT
         % Compute the number of breakpoints, and the minimal and maximal
         % breakpoints of the projection of x + w on the box [xl,xu].
         
         inc = x < xu & w > 0;     % upper bound intersections
         dec = x > xl & w < 0;     % lower bound intersections

         nbrpt = sum(inc | dec);   % number of breakpoints

         % Quick exit if no breakpoints
         if nbrpt == 0
            brptmin = 0;
            brptmax = 0;
            return
         end

         brpt_inc = (xu(inc) - x(inc)) ./ w(inc);
         brpt_dec = (xl(dec) - x(dec)) ./ w(dec);

         brptmin =  inf;
         brptmax = -inf;
         if any(brpt_inc)
            brptmin_inc = min(brpt_inc);
            brptmin = min(brptmin, brptmin_inc);

            brptmax_inc = max(brpt_inc);
            brptmax = max(brptmax, brptmax_inc);
         end
         if any(brpt_dec)
            brptmin_dec = min(brpt_dec);
            brptmin = min(brptmin, brptmin_dec);

            brptmax_dec = max(brpt_dec);
            brptmax = max(brptmax, brptmax_dec);
         end
         
      end % function breakpt
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function sigma = trqsol(x, p, delta)
         %TRQSOL  Largest solution of the TR equation.
         %     This subroutine computes the largest (non-negative) solution
         %     of the quadratic trust region equation
         %
         %           ||x + sigma*p|| = delta.
         %
         %     The code is only guaranteed to produce a non-negative solution
         %     if ||x|| <= delta, and p != 0. If the trust region equation has
         %     no solution, sigma = 0.
         ptx = p'*x;
         ptp = p'*p;
         xtx = x'*x;
         dsq = delta^2;
         
         % Guard against abnormal cases.
         
         rad = ptx^2 + ptp*(dsq - xtx);
         rad = sqrt(max(rad,0));
         
         if ptx > 0
            sigma = (dsq - xtx)/(ptx + rad);
         elseif rad > 0
            sigma = (rad - ptx)/ptp;
         else
            sigma = 0;
         end
         
      end % function trqsol
      
      function [self, flag] = post_iteration(self, x, cgits, successful) %#ok<INUSD>
         % User callback. This routine does nothing by default, but may be
         % overloaded by the user to implement specific functionality.
         flag = 0;
      end
      
   end % methods(Static)
   
end % classdef
