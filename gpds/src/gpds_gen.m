function [x, l] = gpds_gen(base, x, l, gphyp, offset, jit)

  max_rej = 3000;
  batch_sz = 10;
  acc      = 0;
  x_new    = [];
  l_new    = [];

  size_x = length(x);
  it = 0;

  while(~acc)
    x_new     = base(batch_sz);

    if(length(x) > 0)
      K_x       = covSEiso(gphyp, x);
      K_x(1:size(K_x)+1:end) = K_x(1:size(K_x)+1:end) + jit;
      K_x_sqrt  = chol(K_x);
      inv_K_x_sqrt = inv(K_x_sqrt);


      ww = l'*inv_K_x_sqrt;    % Whiten (Zero mean)

      Ks  = covSEiso(gphyp, x, x_new);
      kss = covSEiso(gphyp, x_new);

      pred     = Ks'*inv_K_x_sqrt;
      mn_pred  = pred*ww';
      K_pred   = kss - pred*pred';
  
    else
      mn_pred = 0;
      K_pred  = covSEiso(gphyp, x_new);
    end;
    K_pred(1:size(K_pred)+1:end) = K_pred(1:size(K_pred)+1:end) + jit;

    [K_pred_sqrt, err] = chol(K_pred);
  
    w_new   = randn(size(K_pred_sqrt,1),1);
    l_new   = K_pred_sqrt' * w_new + mn_pred;
    p_new   = sigmoid(l_new,offset);

    for i = 1:batch_sz
      it = it + 1;
      if(rand < p_new(i) | it > max_rej)
        acc = 1;
        break;
      end;
    end;

    x = [x; x_new(1:i)];
    l = [l; l_new(1:i)];
  end;

  x = x(size_x+1:end);
  l = l(size_x+1:end);
