function p_vec = get_p_vec(l_H, G, H, Tmin, rate_scale, gam_sh, gam_sc, offset)

  lmbd = log(sigmoid(l_H,offset));

  G_len = length(G);
  H_len = length(H);

  lg_rs = log(rate_scale);

  prev_event  = Tmin;
  g           = 1;
  next_event  = G(g);

  p_vec = zeros(H_len,1);

  h  = 1;
  while h <= H_len | g <= G_len
    if h > H_len | H(h) > next_event

      if(g <= G_len)
        prev_event      = next_event;

        g = g + 1;
        if(g <= G_len)
          next_event      = G(g);
        else
          next_event      = max(H) + 1;
        end;
      end;
    else
      tmp = -lg_rs+lmbd(h)+log(gamhaz(H(h)-prev_event , gam_sh, gam_sc));
%      p_vec(h) = logdiffexp_v(0,tmp);
      p_vec(h) = tmp;
      h = h + 1;
    end;
  end;
