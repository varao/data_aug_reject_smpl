function [x,l_x,y,l_y] = wr_gpds_gen()

  % What is the prior on the GP mean?
  gp_mean = 0;
  gpm_var = 10;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Log-normal prior on hyperparameters
  gphyp_mn     = log([1;1]);
  gphyp_cvchol = eye(2);

  length_sc_pr = @(theta) gaussian(theta, gphyp_mn, gphyp_cvchol); 

  % Initialize to mean
  gphyp     = gphyp_mn;

  x = [];
  l = [];
  y = [];
  l_y = [];

  for i = 1:10
    [x_new, l_new] = gpds_gen(@(num) randn(num, 1), [x;y], [l;l_y], gphyp);

    y   = [x;x_new(1:end-1)];
    l_y = [l_y;l_new(1:end-1)];
    x   = [x;x_new(end)];
    l_x = [l_x;l_new(end)];
  end;
