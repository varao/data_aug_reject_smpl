function out = nearest_posdef(in)
  in = 0.5 .* (in + in');
  [V,D] = eig(in);
  d     = diag(D);
  thr   = 0.01;
  d(d<thr) = thr;
  D(1:size(D,1)+1:end) = d;
  out = V*D*V';

