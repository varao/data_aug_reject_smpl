function plot_stats(samples,s_g,arg, varargin)

  figure(1);
  clf;
  hold on;
  first = 1;
  last  = length(samples);

  G        = samples(1).G;
  t        = s_g.t;


  v = [];

  for i = first:last
    hyp  = samples(i).gphyp;
    gpm  = samples(i).gp_mean;
    os   = samples(i).offset;
    gam_sh   = samples(i).gam_sh;
    lm   = samples(i).l_max;
    G    = samples(i).G;
    H    = samples(i).H;
    l_G  = samples(i).l_G + gpm;
    l_H  = samples(i).l_H + gpm;
    l_t(i,:)  = s_g.val(i,:);
    l_Gs = lm.*sigmoid(l_G,0);
    l_Hs = lm.*sigmoid(l_H,0);
%    l_ts(i,:) = lm.*sigmoid(l_t(i,:),0);
    xx   = [l_G;l_H];
    xxs  = [l_Gs;l_Hs];

    if nargin == 2 || arg == 1
      plot(t,l_t(i,:),'r','markersize',.3);
    elseif arg == 2
      v = [v,hyp(1)];
    elseif arg == 3
      v = [v,hyp(2)];
    elseif arg == 4
      v = [v,gpm];
    elseif arg == 5
      v = [v,lm];
    elseif arg == 6
      v = [v,os];
    elseif arg == 7
      v = [v,gam_sh];
    end;
  end

  if nargin == 2 || arg == 1

    plot(G,0.*G, 'b*');
    plot(t,mean(l_t,1),'g');
%   for j = 1:num_bins
%       mn(j) = mean(smpl{j,:});
%       vr(j) = var(smpl{j,:});
%   end;
%   plot(bins,mn,'r');
%   hold on;
%   plot(bins,mn+ sqrt(vr),'r-','markersize',1);
%   plot(bins,mn- sqrt(vr),'r-','markersize',1);
%   plot(G, mean([max(mn+sqrt(vr)),min(mn-sqrt(vr))]) ,'b.');
%    ylim([0,max(mn+sqrt(vr))+.1]);
%   ylim([min(mn-sqrt(vr)),max(mn+sqrt(vr))+.1]);
%   xlim([0,max(G)+5]);
%   ylabel('Intensity');
%   title('Coal mine data (Shape parameter = 2)');
%   set(gcf,'Paperposition',[1 1 8 6])
%   print -depsc tmp.eps
%   !epspdf tmp.eps
  else
    if arg == 2
      display 'Length scale';
    elseif arg == 3
      display 'Variance';
    elseif arg == 4
      display 'GP mean';
    elseif arg == 5
      display 'Lmax';
    elseif arg == 6
      display 'Offset';
    elseif arg == 7
      display 'Gamma';
    end;
    if nargin == 3
      hist(v,20);
    else
      plot(v);
    end;
  end;


