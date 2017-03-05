clf
pt = plot_dens(c,d,pos*40,'-');
pt1 = plot_dens(a,b,neg*40,'--');
%legend([pt,pt1],'positive','control')
ylim([0,.65])
%plot(pos*40,0,'k+','markers',3)
%ylim([-5,10.65])
box on;
print -depsc plot_galaxy.eps
!epstopdf plot_galaxy.eps
