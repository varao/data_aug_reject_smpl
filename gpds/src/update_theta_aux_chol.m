function [theta, ff, U] = update_theta_aux_chol(theta, ff, Lfn, Kfn, theta_Lprior, slice_width, U)
%UPDATE_THETA_AUX_CHOL MCMC update to GP hyperparam. Fixes nu used to draw f, rather than f itself
%
% Lfn @fn Log-likelihood function, Lfn(ff) returns a scalar

% Iain Murray, November 2009

Ufn = @(th) chol(Kfn(th));
DEFAULT('theta_Lprior', @(l) log(double((l>log(0.1)) && (l<log(10)))));
DEFAULT('U', Ufn(theta));

% Change of variables
nu = U' \ ff;

% Slice sample theta|nu
particle = struct('pos', theta);
particle = eval_particle(particle, -Inf, nu, Lfn, theta_Lprior, U);
step_out = (slice_width > 0);
slice_width = abs(slice_width);
slice_fn = @(pp, Lpstar_min) eval_particle(pp, Lpstar_min, nu, Lfn, theta_Lprior, Ufn);
particle = slice_sweep(particle, slice_fn, slice_width, step_out);
theta = particle.pos;
ff = particle.ff;
U = particle.U;


function pp = eval_particle(pp, Lpstar_min, nu, Lfn, theta_Lprior, U)
% U is a precomputed chol(Kfn(pp.pos)) or a function that will compute it

% Prior
theta = pp.pos;
Ltprior = theta_Lprior(theta);
if Ltprior == -Inf
    pp.on_slice = false;
    return;
end

if ~isnumeric(U)
    U = U(theta);
end
ff = (nu'*U)';

pp.Lpstar = Ltprior + Lfn(ff);
pp.on_slice = (pp.Lpstar >= Lpstar_min);
pp.U = U;
pp.ff = ff;
