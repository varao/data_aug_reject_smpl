function ll = log_sig(l, G, H, Tmin, rate_scale, gam_sh, gam_sc, offset)

  % Log-lik under this a ligmoid likelihood function

  lmbd = log(sigmoid(l,offset));

  G_len = length(G);
  H_len = length(H);

  l_G   = lmbd(1:G_len);
  l_H   = lmbd(G_len+1:end);

  lg_rs = log(rate_scale);

  prev_event  = Tmin;
  g           = 1;
  next_event  = G(g);

  lg1 = zeros(G_len,1);
  lg2 = zeros(H_len,1);

  h  = 1;
  while h <= H_len | g <= G_len
    if h > H_len | H(h) > next_event

      if(g <= G_len)
        lg1(g) = (next_event-prev_event);
        prev_event      = next_event;

        g = g + 1;
        if(g <= G_len)
          next_event      = G(g);
        else
          next_event      = max([H;G]) + 1;
        end;
      end;
    else
      lg2(h) = (H(h)-prev_event);
      h      = h + 1;
    end;
  end;

  p2 = sum(l_G) - G_len*lg_rs + sum(log(gamhaz(lg1,gam_sh, gam_sc)));
  if(sum(isinf(p2))) 
    error 'WTF?';
  end;

  p = gamhaz(lg2,gam_sh, gam_sc);
  p = 1 - exp(l_H-lg_rs).*p;

  p = sum(log(p));

  ll = p + p2;
