function med = plot(sample, rest, galaxy)

  Tmin = 0;
  Tmax =  5;
  t = linspace(Tmin, Tmax, 200)';
  num_smpl = length(galaxy);

  t = rest.t(:)';
  for s = 1:length(sample)

    if(s == 1)
      G_hist = sample(s).G(num_smpl:end);
    else
      G_hist = [G_hist;sample(s).G(num_smpl:end)];
    end;
    len(s) = length(sample(s).G)-num_smpl;
  end;
  clf;
  size(G_hist)
%  gr = [-2:.01:12];
%  k = ksdensity(G_hist, gr );
%  plot(gr,k,'linewidth',5);
%  hold on;
%  plot(galaxy(1:num_smpl),0,'+k');
   ed = [0:30:600];
   hc = histc(len,ed);
   max(len);
   med = bar(ed,hc./sum(hc));
   xlim([-30,600]);
   ylabel('Probability');
   xlabel('Number of rejections');
%  xlabel('Normalized velocity');
%  ylabel('Density of rejections');
  set(gcf,'Paperposition',[1 1 3 2]);
  print -depsc plot_thin2.eps
  !epstopdf plot_thin2.eps
