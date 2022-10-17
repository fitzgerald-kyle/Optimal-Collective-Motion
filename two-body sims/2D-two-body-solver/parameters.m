function params = parameters(nmax)

% This function contains parameters that are used when solving the BVP

params.nmax = nmax; % maximum number of iterations to run BVP loop
params.tol = 1e-6; % error tolerance for solution of BVP
params.maxstep = 1; % maximum step size to make in p0 when solving BVP
params.step_ratio = 0.5; % step size decrease to use in line search
params.decrease_param = 0; % minimum error decrease in line search
params.minstep = 1e-6; % minimum step size in line search
params.ode_options = odeset('AbsTol',1e-11,'RelTol',1e-11,'Events',@conj_pt_test); % parameters for solving ODEs

end