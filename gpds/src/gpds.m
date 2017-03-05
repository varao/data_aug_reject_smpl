function [sample, vargout] = gpds(data, resln, varargin)

  if nargin == 2

    grid = 1;
    x    = data;

    burnin      = 500;
    num_samples = 2000;
    aux = .1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hyperparameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    offset  = 0;

    % What is the prior on the GP mean?
    gp_mean = 0;
    gpm_var = 20;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log-normal prior on hyperparameters
    gphyp_mn     = log([1;10]);
    gphyp_mn(2)  = gphyp_mn(2)/2;
    gphyp_cvchol = 1.* eye(2);
    gphyp_cvchol(1,1) = 1;

    length_sc_pr = @(theta) gaussian(theta, gphyp_mn, gphyp_cvchol); 

    % Initialize to mean
    gphyp     = log([1;2]);
%    gphyp     = gphyp_mn;

    Tmin =  0;
    Tmax =  5;
    if grid == 1
      t = linspace(Tmin, Tmax, 200)';
      sample_grid.t = t;
    end;
  else

    x    = data.eventtimes;
    Tmin = data.Tmin;
    Tmax = data.Tmax;

    mcmc_pars = varargin{1};

    burnin      = mcmc_pars.burnin;
    num_samples = mcmc_pars.num_samples;

    aux         = 1;

    grid = mcmc_pars.grid;
    if grid == 1
      t = linspace(Tmin, Tmax, 20)';
      sample_grid.t = t;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hyperparameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hyper_pars  = varargin{2};

    offset  = 0;

    % What is the prior on the GP mean?
    gp_mean = hyper_pars.gp_mean;
    gpm_var = hyper_pars.gp_var;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log-normal prior on hyperparameters
    gphyp_mn     = hyper_pars.gphyp_mn;
    gphyp_cvchol = hyper_pars.gphyp_cvchol;

    % Actually, should just pass a function handle
    length_sc_pr = @(theta) gaussian(theta, gphyp_mn, gphyp_cvchol); 

    % Initialize to mean
    gphyp     = gphyp_mn;

  end;
  offset  = 0;

  clf;
  plot(x,0,'*b');
  hold on;

  base_mean = mean(data);
  base_std  = sqrt(var(data));
  base = @(num) base_std*randn(num, 1) + base_mean;
%  base = @(num) rand(num,1).*10;
  G   = x;
  num_obs = length(x);

  mat_d = 5;
  jit = 1e-2;
  log_jit = log(jit);
%  K_G = covMaterniso_jit(mat_d, gphyp, G, jit);
  K_G = covSEiso_jit(gphyp, G, jit);
  K_G_sqrt = chol(K_G);  % A'*A = B, Chol returns A, an upper triangular matrix

  w_G = randn(size(K_G_sqrt,1),1);
  l_G = K_G_sqrt' * w_G;   % Has zero mean
  l_x = l_G;

  for(s = -burnin:num_samples)

    y_new = []; l_y_new = [];
    for r = 1:num_obs
      [y_tmp, l_tmp] = gpds_gen(base, G, l_G, gphyp, offset, jit);

      if(length(y_tmp) > 1)
        y_new   = [y_new; y_tmp(1:end-1)];
        l_y_new = [l_y_new; l_tmp(1:end-1)];
      end;
    end;

    l_x = l_G(1:num_obs);
    y   = y_new;
    l_y = l_y_new;
    G   = [x;y];
    l_G = [l_x;l_y];
    ff = l_G;

    if(rem(s,10) == 0)
      s
      gphyp'
      [offset base_mean, base_std]
      length(l_y)
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GP hyperparameters
    % Assuming l_G, l_H and ff have 0-mean

%   [gphyp, ff, aux, K_G_sqrt] = update_theta_aux_surr(gphyp, l_G, @(x) gpds_lik(x, num_obs, offset), ...
%                                                               @(x) covMaterniso_jit(mat_d, x, G,jit), aux, length_sc_pr, 1);
%   if(rem(s,10) == 0)
%     [gphyp, ff, aux, K_G_sqrt] = update_theta_aux_surr(gphyp, l_G, @(x) gpds_lik(x, num_obs, offset), ...
%                                                               @(x) covSEiso_jit(x, G,jit), aux, length_sc_pr, 1);
%   end;                                                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    K_G = covMaterniso_jit(mat_d, gphyp, G, jit);
    K_G = covSEiso_jit(gphyp, G, jit);
    K_G_sqrt = chol(K_G);  % A'*A = B, Chol returns A, an upper triangular matrix
    inv_K_G_sqrt = inv(K_G_sqrt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot grid
    if s>0 & grid == 1
%     Kt  = covMaterniso(mat_d, gphyp, G, t);
%     ktt = covMaterniso(mat_d, gphyp, t);
      Kt  = covSEiso(gphyp, G, t);
      ktt = covSEiso(gphyp, t);

      pred_t     = Kt'*inv_K_G_sqrt;
      mn_pred_t  = pred_t*inv_K_G_sqrt'*ff;
      K_pred_t   = ktt - pred_t*pred_t';
      K_pred_t(1:size(K_pred_t)+1:end) = K_pred_t(1:size(K_pred_t)+1:end) + jit;
  
      K_pred_t_sqrt = chol(nearest_posdef(K_pred_t));
      l_t = mn_pred_t + K_pred_t_sqrt' * randn(size(K_pred_t_sqrt,1),1);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curr_log_like = gpds_lik([l_x; l_y], num_obs, offset);

    w_pr = randn(size(K_G_sqrt,1),1);
    l_pr = K_G_sqrt' * w_pr;

    if(s == 200)
      s
    end;
  
    [l_G, phi, curr_log_like] = elliptical_slice(l_G, l_pr, @(x) gpds_lik(x, num_obs, offset), curr_log_like, 0);

    dist_offset = @(offset) gpds_lik(x, num_obs, offset) - (offset^2)/2;
    offset = slicesample(offset ,1,'logpdf', dist_offset);

    smpl    = [x;y];
    num_smp = length(smpl);
    sum_smp = sum(smpl);
    avg_smp = sum_smp/num_smp;
    sum_sq_smp = sum((smpl-avg_smp).^2);

    nu = 0.1; 

    nig_p1  = sum_smp / (nu + num_smp);
    nig_p2  = nu + num_smp;
    nig_p3  = 1 + num_smp/2;
    nig_p4  = 10 + 0.5 * sum_sq_smp + 0.5 * num_smp * (avg_smp^2) / (num_smp/nu + 1);
   
    base_pr   = gamrnd(nig_p3, 1/nig_p4);
    base_var  = 1/base_pr;
    base_mean = sqrt(base_var/nig_p2) * randn + nig_p1;
    base_std = sqrt(base_var);

    base = @(num) base_std*randn(num, 1) + base_mean;

%   hold on;
    if(s > 0)
%     subplot(2,1,1);
%     plot(t,(l_t + gp_mean),'m');
%     subplot(2,1,2);
%     plot(t,l_max.*sigmoid(l_t+gp_mean,offset),'m');


      sample(s).G = G;
      sample(s).l_G = l_G;
      sample(s).gphyp = gphyp;
      sample(s).gp_mean = gp_mean;
      sample(s).offset = offset;
      sample(s).base_mean = base_mean;
      sample(s).base_std = base_std;

      if grid == 1
        if(s == 1)
          s_mean = sigmoid(l_t,offset);
        else
          s_mean = s_mean + sigmoid(l_t,offset);
        end;
        sample_grid.val(s,:) = l_t;
      end;
      plot(t, normpdf(t, base_mean, base_std) .* sigmoid(l_t, offset) , 'g');
    end;
  end;
  plot(t, normpdf(t, base_mean, base_std) .* s_mean./s);
  plot(t, s_mean./s, 'g');
  hold on;
  plot(t, normpdf(t, base_mean, base_std) .* mean(s_mean./s), 'r');

  vargout(1) = sample_grid;

function K_GH = covMaterniso_jit(mat_d, gphyp, GH,jit)
    K_GH = covMaterniso(mat_d, gphyp, GH);
    K_GH(1:size(K_GH)+1:end) = K_GH(1:size(K_GH)+1:end) + jit;

function K_GH = covSEiso_jit(gphyp, GH,jit)
    K_GH = covSEiso(gphyp, GH);
    K_GH(1:size(K_GH)+1:end) = K_GH(1:size(K_GH)+1:end) + jit;

function varargout = gaussian(x, mn, cholSigma)

    % Return the log probability of x.
    varargout(1) = {mvnormpdfln(x, mn, cholSigma)};

    % Perhaps also return the derivative.
    if nargout > 1
      varargout(2) = {-solve_triu(cholSigma, ...
                                  solve_tril(cholSigma',x-mn))};

    end

