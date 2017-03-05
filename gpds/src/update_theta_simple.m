function [theta, U] = update_theta_simple(theta, ff, Lfn, Kfn, theta_Lprior, slice_width, U)
%UPDATE_THETA_SIMPLE Standard slice-sampling MCMC update to GP hyper-param

% Iain Murray, January 2010


Ufn = @(th) chol(Kfn(th));
DEFAULT('theta_Lprior', @(l) log(double((l>log(0.1)) && (l<log(10)))));
DEFAULT('U', Ufn(theta));

% Slice sample theta|ff
particle = struct('pos', theta, 'ff', ff);
particle = eval_particle(particle, -Inf, Lfn, theta_Lprior, U);
step_out = (slice_width > 0);
slice_width = abs(slice_width);
slice_fn = @(pp, Lpstar_min) eval_particle(pp, Lpstar_min, Lfn, theta_Lprior, Ufn);
particle = slice_sweep(particle, slice_fn, slice_width, step_out);
theta = particle.pos;
U = particle.U;

function pp = eval_particle(pp, Lpstar_min, Lfn, theta_Lprior, U)

% Prior
theta = pp.pos;
Ltprior = theta_Lprior(theta);

if Ltprior == -Inf
    pp.Lpstar   = -Inf;
    pp.on_slice = false;
    return;
end

if ~isnumeric(U)
    U = U(theta);
end

Lfprior = -0.5*(pp.ff'*solve_chol(U, pp.ff)) - sum(log(diag(U))); % + const
pp.Lpstar = Ltprior + Lfprior + Lfn(pp.ff);
pp.on_slice = (pp.Lpstar >= Lpstar_min);
pp.U = U;
