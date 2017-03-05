function G = mod_renew_func(rate, rate_max, Tmin, Tmax, gam_sh)

  gam_sc  = 1/gam_sh;
  haz_max = gam_sh;

  Omega = rate_max * haz_max;

  M = poissrnd(Omega*(Tmax-Tmin));
  T = rand(M,1).*(Tmax-Tmin) + Tmin;
  T = sort(T);
  T = [Tmin;T];

  lmbd = rate(T)./rate_max;

  last_event = 1;
  t = linspace(Tmin, Tmax, 200)';

  G = [];

  for i = 1:M
    p = lmbd(i)*gamhaz(T(i+1)-T(last_event) , gam_sh, gam_sc)/haz_max;
    U=rand;
    if U <= p
      last_event = i+1;
      G = [G;T(i+1)];
    end
  end

  T = T(2:end);
  l = rate(t);

  clf;
  plot(t,l)
  hold on;
  plot(G,G.*0,'*r')

end
