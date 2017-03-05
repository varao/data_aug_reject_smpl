function [sample, vargout] = mod_renew_inf(data, varargin)

  if nargin == 1

    G    = data;
    grid = 1;

    burnin      = 500;
    num_samples = 2000;
    aux = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hyperparameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%
    % Does the sigmoid have an offset?  Actually, use a nonzero GP mean instead
     offset  = 2;
    %offset_var = 10;
    %offset_prop_std  = 1;

    % What is the prior on the GP mean?
    gp_mean = 0;
    gpm_var = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Log-normal prior on hyperparameters
    gphyp_mn     = log([1;1]);
    gphyp_cvchol = eye(2);

  %  length_sc_pr = @(theta) -sum((theta-gphyp_mn).^2)/2; 
    length_sc_pr = @(theta) gaussian(theta, gphyp_mn, gphyp_cvchol); 

    % Initialize to mean
    gphyp     = gphyp_mn;

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Renewal time distribution

    gam_sh_prop_std  = .1;
    gam_sh  = 1;
    gam_sc  = 1/gam_sh;
    haz_max = gam_sh;

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Prior on scaling of sigmoid ie on l_max
    l_sh = 5;
    l_sc = 1;

    Tmin = 0;
    Tmax = max(G);
    if grid == 1
      t = linspace(Tmin, Tmax, 200)';
      sample_grid.t = t;
    end;
  else

    G    = data.eventtimes;
    Tmin = data.Tmin;
    Tmax = data.Tmax;

    mcmc_pars = varargin{1};

    burnin      = mcmc_pars.burnin;
    num_samples = mcmc_pars.num_samples;

    aux         = 1;

    grid = mcmc_pars.grid;
    if grid == 1
      t = linspace(Tmin, Tmax, 200)';
      sample_grid.t = t;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Hyperparameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hyper_pars  = varargin{2};

    % Does the sigmoid have an offset?  Actually, use a nonzero GP mean instead
     offset  = 2;

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

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Renewal time distribution

    gam_sh_prop_std  = hyper_pars.gam_sh_prop_std;
    gam_sh           = hyper_pars.gam_sh;
    gam_sc           = 1/gam_sh;
    haz_max          = gam_sh;

    %%%%%%%%%%%%%%%%%%%%%%%%
    % Prior on scaling of sigmoid ie on l_max
    l_sh = hyper_pars.l_sh;
    l_sc = hyper_pars.l_sc;
  end;

  len_G = length(G);
  figure(1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  jit = 1e-5;
  log_jit = log(jit);

  K_G = covSEiso_jit(gphyp, G, jit);
  K_GH_sqrt = chol(K_G);  % A'*A = B, Chol returns A, an upper triangular matrix
  w_G = randn(size(K_GH_sqrt,1),1);
  l_G = K_GH_sqrt' * w_G;   % Has zero mean
  
% clf;
% subplot(2,1,1);
% plot(G,0,'b*');
% hold on;
% subplot(2,1,2);
% plot(G,0,'b*');
% hold on;

  H   = [];
  w_H = [];
  l_H = [];

  l_max = gamrnd(l_sh,l_sc);

  for(s = -burnin:num_samples)

    Omega = l_max * haz_max;
    Rate = Omega*(Tmax-Tmin);

    if(rem(s,20) == 0)
      s
      gam_sh
      l_max
      gphyp'
      gp_mean
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GP hyperparameters
    % Assuming l_G, l_H and ff have 0-mean
    if(rem(s,1) == 1)
      [gphyp, ff, aux, K_GH_sqrt] = update_theta_aux_surr(gphyp, [l_G;l_H], @(x) log_sig(x+gp_mean, G, H, Tmin, haz_max, gam_sh, gam_sc, offset), ...
                                                                @(x) covSEiso_jit(x, [G;H],jit), aux, length_sc_pr, 1);
    else
      ff = [l_G;l_H];
      K_GH_sqrt = chol(covSEiso_jit(gphyp, [G;H], jit));
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%
    %   [gphyp, K_GH_sqrt] = update_gphp_hmc(gphyp, [G;H], [l_G;l_H], {'covSum', {'covSEiso', 'covNoise'}}, log_jit , length_sc_pr, 10, 0.001);
    %   ff = [l_G;l_H];

    %%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    inv_K_GH_sqrt = inv(K_GH_sqrt);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot grid
    if s>0 & grid == 1
      Kt  = covSEiso(gphyp, [G;H], t);
      ktt = covSEiso(gphyp, t);

      pred_t     = Kt'*inv_K_GH_sqrt;
      mn_pred_t  = pred_t*inv_K_GH_sqrt'*ff;
      K_pred_t   = ktt - pred_t*pred_t';
      K_pred_t(1:size(K_pred_t)+1:end) = K_pred_t(1:size(K_pred_t)+1:end) + jit;
  
      K_pred_t_sqrt = chol(K_pred_t);
      l_t = mn_pred_t + K_pred_t_sqrt' * randn(size(K_pred_t_sqrt,1),1);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Resample mean

%    l_G = ff(1:len_G) + gp_mean;
%    l_H = ff(len_G+1:end) + gp_mean;


%    ww = ([l_G;l_H])'*inv_K_GH_sqrt;    % ww' actually, and doesn't have 0 mean
%    I_isq_K_GH = sum(inv_K_GH_sqrt,1);

%   mn_pr = I_isq_K_GH*I_isq_K_GH' + 1/gpm_var;
%   mn_mn = (ww*I_isq_K_GH')/mn_pr;

%    gp_mean = mn_mn + randn/sqrt(mn_pr);
%   [mn_mn,1/mn_pr, gphyp]
%    curr_log_like = log_sig([l_G;l_H]+ gp_mean, G, H, Tmin, haz_max, gam_sh, gam_sc, offset);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Block Gibbs resample thinned events
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    M     = poissrnd(Rate);

    if M > 0
      H_new = rand(M,1).*(Tmax-Tmin) + Tmin;
      H_new = sort(H_new);

      ww = ff'*inv_K_GH_sqrt;    % Whiten (Zero mean)

      Ks  = covSEiso(gphyp, [G;H], H_new);
      kss = covSEiso(gphyp, H_new);

      pred     = Ks'*inv_K_GH_sqrt;
      mn_pred  = pred*ww';
      K_pred   = kss - pred*pred';
      K_pred(1:size(K_pred)+1:end) = K_pred(1:size(K_pred)+1:end) + jit;
  
      [K_pred_sqrt, err] = chol(K_pred);
  
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Overwrite old stuff

      l_G = ff(1:len_G);

      w_H   = randn(size(K_pred_sqrt,1),1);
      l_H   = K_pred_sqrt' * w_H + mn_pred;

      p_vec =  get_p_vec(l_H+gp_mean, G, H_new, Tmin, Omega/l_max, gam_sh, gam_sc, offset);
      u_H   =  rand(length(p_vec),1);
      thin  =  p_vec < log(u_H);
      H     =  H_new(thin);
      w_H   =  w_H(thin);
      l_H   =  l_H(thin);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end;

    curr_log_like = log_sig([l_G;l_H]+ gp_mean, G, H, Tmin, haz_max, gam_sh, gam_sc, offset);

    w_pr1 = randn(size(K_GH_sqrt(1:len_G,1:len_G),1),1);
    l_pr1 = K_GH_sqrt(1:len_G,1:len_G)' * w_pr1;
  
    if M > 0
      inv_K_G_sqrt = inv(K_GH_sqrt(1:len_G,1:len_G));
      pred         = Ks(1:len_G,:)'*inv_K_G_sqrt;
      w_pr2        = randn(size(K_pred_sqrt,1),1);
      l_pr2        = K_pred_sqrt' * w_pr2 + pred*w_pr1;
      w_pr2 = w_pr2(thin);
      l_pr2 = l_pr2(thin);
    else
      l_pr2 = [];
    end;

    [xx, phi, curr_log_like] = elliptical_slice([l_G;l_H], [l_pr1;l_pr2], @(x) log_sig(x+gp_mean, G, H, Tmin, haz_max, gam_sh, gam_sc, offset), curr_log_like, 0);

    l_G = xx(1:len_G);
    l_H = xx(len_G+1:end);

%   offset_new = offset + offset_prop_std*randn();
%   ll_new     = log_sig([l_G;l_H], G, H, Tmin, haz_max, gam_sh, gam_sc, offset_new);
%   acc        = (ll_new - (offset_new^2)/(2*offset_var)) - (curr_log_like - (offset^2)/(2*offset_var));

%  if(log(rand) < acc)
%    offset        = offset_new;
%    curr_log_like = ll_new;
%  end;

    for gl = 1:10
      gam_sh_new = gam_sh + gam_sh_prop_std*randn();

      if(gam_sh_new >= 1)
        ll_old     = log_sig([l_G;l_H], G, H, Tmin, gam_sh, gam_sh, 1/gam_sh, offset);
        ll_new     = log_sig([l_G;l_H], G, H, Tmin, gam_sh_new, gam_sh_new, 1/gam_sh_new, offset);

        if(~isreal(ll_new))
          error 'WTF?'
        end;
        acc        = (ll_new - (gam_sh_new-1)/10 - (ll_old - (gam_sh-1)/10));

%      [curr_log_like ll_old ll_new acc]
%      [gam_sh gam_sh_new]

        if(log(rand) < acc)
          gam_sh        = gam_sh_new;
          gam_sc  = 1/gam_sh;
          haz_max = gam_sh;
          curr_log_like = ll_new;
%         display('Accept!');
        end;
      end;
    end;

    inv_l_sc_post = Rate/l_max + 1/l_sc;
    l_max = gamrnd(l_sh+length(xx), 1/inv_l_sc_post);
    

    if(s > 0)
%     [tt, tt_i] = sort([G;H]);
%     subplot(2,1,1);
%%      plot(tt,(xx(tt_i)+gp_mean),'m');
%     plot(t,(l_t + gp_mean),'m');
%     subplot(2,1,2);
%%      plot(tt,l_max.*sigmoid(xx(tt_i)+gp_mean,offset),'m');
%     plot(t,l_max.*sigmoid(l_t+gp_mean,offset),'m');


      sample(s).H = H;
      sample(s).G = G;
      sample(s).l_H = l_H;
      sample(s).l_G = l_G;
      sample(s).l_max = l_max;
      sample(s).gphyp = gphyp;
      sample(s).gp_mean = gp_mean;
      sample(s).offset = offset;
      sample(s).gam_sh = gam_sh;

      if grid == 1
        sample_grid.val(s,:) = l_max .* sigmoid(l_t,offset);
      end;
    end;
  end;

  vargout(1) = sample_grid;

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


