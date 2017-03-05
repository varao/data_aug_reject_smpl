tlim_arr = [200 100 50 20 10 5];

for i = 1:length(tlim_arr)
  G_temp = G(G< tlim_arr(i));
  size(G_temp)
  [sample, misc] = mod_renew_inf(G_temp');
  rslt{i} = sample;
end;
