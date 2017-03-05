function [theta, cholK] = update_gphp_hmc(theta, xx, ff, Kfn, jit, gp_hyp_prior, tau, epsilon)

debug_flg = 0;

% Calculate the negative log marginal likelihood and gradient.
[nlml dnlml] = gpr([theta;jit], Kfn, xx, ff);

% Calculate the prior and gradient.
[logprior dlogprior] = gp_hyp_prior(theta);

% Instantiate momenta.
P = randn(size(theta));

% Calculate the Hamiltonian.
init_H = P'*P/2 + nlml - logprior;

% Calculate the gradient.
G = dnlml(1:end-1) - dlogprior;

% Initialize
loghp = theta;

% Take tau steps.
for t=1:tau

  % Take a half-step in momentum space.
  P = P - epsilon * G/2;

  % Take a step in parameter space.
  loghp = loghp + epsilon * P;

  if debug_flg == 1
    % Calculate the distance.
    distance = sqrt(sum((loghp-theta).^2));
    fprintf('%d - Distance %f\n', t, distance);
  end

  % Reevaluate the gradient.
  [nlml dnlml]         = gpr([loghp; jit], ...
                             Kfn, xx, ff);
  [logprior dlogprior] = gp_hyp_prior(loghp);
  G = dnlml(1:end-1) - dlogprior;

  % Take another half-step in momentum space.
  P = P - epsilon * G/2;
end

% Evaluate the new Hamiltonian.
[nlml dnlml]         = gpr([loghp; jit], Kfn, xx, ff);
[logprior dlogprior] = gp_hyp_prior(loghp);

new_H = P'*P/2 + nlml - logprior;

if debug_flg == 1
  fprintf('Old Ham: %f  New Ham: %f  Prob: %f\n', init_H, new_H, ...
          exp(init_H - new_H));
end

% Accept or reject.
if rand() < exp(init_H - new_H)

  % Update the state.
  theta = loghp;
end

  % Recompute the Cholesky decomposition.
  cholK = chol(feval(Kfn{:}, ...
                        [theta; jit] , xx));

end % update_gphp



