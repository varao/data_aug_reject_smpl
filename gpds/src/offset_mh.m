   offset_new = offset + offset_prop_std*randn();
   ll_new     = log_sig([l_G;l_H], G, H, Tmin, Omega/l_max, gam_sh, gam_sc, offset_new);
   acc        = (ll_new - (offset_new^2)/(2*offset_var)) - (curr_log_like - (offset^2)/(2*offset_var));

   if(log(rand) < acc)
     offset        = offset_new;
     curr_log_like = ll_new;
   end;

