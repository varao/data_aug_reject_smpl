function p = renewal_prob(G, t, l_t, hazard, varargin);

  if nargin == 5
    l_t = varargin{1} .* sigmoid(l_t,0);
  end;

  l_G = log(interp1(t, l_t, G));

  wait_time = [G;t(end)] - [t(1);G];

  tau   = zeros(1,length(t));

  tprev = t(1);
  for g = [G;t(end)]'
    indx       = (t > tprev) & (t <= g);
    tau(indx)  = t(indx)-tprev;
    if min(tau) < 0
      error 'WTF?';
    end;
    tprev      = g;
  end;

  int1 = trapz(t, l_t.* hazard(tau));
  int1 = sum(l_G) + sum(log(hazard(wait_time))) - int1; 
  
  int1 = int1 - log(l_t(end)) - log(hazard(wait_time(end)));

  p = int1;
