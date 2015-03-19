
function opts = numerics_options(opts)

% nrsolve options
opts.nr.max_iterations = 20;
opts.nr.tolerance = 1e-9;
opts.nr.alpha_min = 0;          % the minimum step size for the newton step
opts.nr.mu = 0;                 % the parameter for the armijo condition
opts.nr.verbose = 0;
opts.nr.linesearch = 'backtrack';

% time domain integration options
opts.sim.max_iters = 40;
opts.sim.tolerance = 1e-9; 
opts.sim.var_step = true;
opts.sim.eps_thresh = 10^-4;

% the following invalidates all of the above
opts.nr.use_fsolve = false;
