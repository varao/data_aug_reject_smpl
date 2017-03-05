function med = plot_dens(sample, rest, galaxy, lnsty)

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
      s_mean = sigmoid(l_t,offset);
    else
      s_mean = s_mean + sigmoid(l_t,offset);
    end;
    pdf(s,:) = normpdf(t, base_mean, base_std) .* sigmoid(l_t, offset);
    pdf(s,:) = 20 .* pdf(s,:) ./ sum(pdf(s,:));
%    pdf(s,:) = l_t-offset;
  end;
 % clf;
 %plot(t, median(pdf,1));
  hold on;
  bot = prctile(pdf,10);
  top = prctile(pdf,90);
  X = [t, fliplr(t)];
  Y = [bot, fliplr(top)];
% h = fill(X,Y,'b');
% set(h,'facealpha',.1);
% plot(t, bot,'b', 'linewidth',2);
% plot(t, top,'b', 'linewidth',2);
% plot(t, median(pdf,1), 'linewidth',2,'color','r');
  plot(t, bot, lnsty, 'linewidth',1,'color',[0 0 0]+.5);
  plot(t, top, lnsty, 'linewidth',1,'color',[0 0 0]+.5);
  pt = plot(t, median(pdf,1), lnsty, 'linewidth',1,'color','k');
%  plot(galaxy,0,'+k')
  xlabel('Normalized voice shimmer');
  ylabel('Probability density');
  ylim([-.00,0.9 .* max(max(pdf))]);
%  legend('positive','negative')
  set(gca,'xlim',[Tmin,Tmax]);
  set(gca,'xtick',[Tmin:Tmax]);
  set(gcf,'Paperposition',[1 1 3 2]);
  colormap(gray);
  print -depsc plot_galaxy.eps
  !epstopdf plot_galaxy.eps
  med = median(pdf,1);
  med(85:95);
  t(85:95);
  med = pt;

