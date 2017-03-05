function [G,l,T] = mod_renew_prior()

  gam_sh  = 2;
  gam_sc  = 1/gam_sh;
  haz_max = gam_sh;

  l_max   = 5;

  Omega = l_max * haz_max;

  Tmin = 0;
  Tmax = 200;

  M = poissrnd(Omega*(Tmax-Tmin));
  T = rand(M,1).*(Tmax-Tmin) + Tmin;
  T = sort(T);

  K = covSEiso(log([1,1]), T);
  K(1:size(K,1)+1:end) = K(1:size(K,1)+1:end) + 1e-5;
  T = [Tmin;T];

% l = mvnrnd(zeros(1,M),K);
  [K_sqrt,err] = chol(K);
  if err ~= 0
      error('Bad covariance matrix.');
  end
  l = K_sqrt' * randn(size(K_sqrt,1),1);

  l = l - mean(l) + 1;

  lmbd = sigmoid(l,2);

  last_event = 1;

  G = [];

  for i = 1:M
    p = lmbd(i)*gamhaz(T(i+1)-T(last_event) , gam_sh, gam_sc)/haz_max;
    U=rand;
    if U <= p
      last_event = i+1;
      G = [G,T(i+1)];
    end
  end

  T = T(2:end);

  l = l_max.*lmbd;
  clf;
  plot(T,l)
%  plot(T,l_max.*lmbd,'g')
  hold on;
  plot(G,G.*0,'*r')
end
