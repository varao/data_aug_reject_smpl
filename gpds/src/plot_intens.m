function med = plot_intens(sample, rest, galaxy)

  Tmin = -0;
  Tmax =  5;
  t = linspace(Tmin, Tmax, 200)';

  t = rest.t(:)';
  for s = 1:length(sample)
    l_t = rest.val(s,:);
    gphyp = sample(s).gphyp;
    gp_mean = sample(s).gp_mean;
    offset = sample(s).offset;
    base_mean = sample(s).base_mean;
    base_std = sample(s).base_std;

    if(s == 1)
      s_mean = (l_t-offset);
    else
      s_mean = s_mean + (l_t-offset);
    end;
    pdf(s,:) = (l_t- offset);
  end;
  clf;
  plot(t, median(pdf,1));
  hold on;
  bot = prctile(pdf,10);
  top = prctile(pdf,90);
  X = [t, fliplr(t)];
  Y = [bot, fliplr(top)];
%  h = fill(X,Y,'b');
% set(h,'facealpha',.1);
% plot(t, bot,'b', 'linewidth',2,'linestyle','-');
% plot(t, top,'b', 'linewidth',2,'linestyle','-');
% plot(t, median(pdf,1), 'linewidth',2,'color','r');
  plot(t, bot, 'linewidth',2,'color',[0 0 0]+.5);
  plot(t, top, 'linewidth',2,'color',[0 0 0]+.5);
  plot(t, median(pdf,1), 'linewidth',2,'color','k');
  plot(galaxy,0,'+k')
  xlabel('Normalized voice shimmer');
  ylabel('Modulating intensity');
  ylim([-6,10]);
  set(gca,'xlim',[Tmin,Tmax]);
  set(gca,'xtick',[Tmin:Tmax]);
  set(gcf,'Paperposition',[1 1 3 2])
  print -depsc plot_galaxy_int.eps
  !epstopdf plot_galaxy_int.eps
  med = median(pdf,1);
  med(85:95);
  t(85:95);
