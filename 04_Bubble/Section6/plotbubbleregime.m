[lik, probsT, At_mat, ~, bubble] = evallp_mod(phi,dds,gx,hx,nstates,ETAMATRIX,zzd) ;

% time = 1994:.00:2016.25 ;
time = 1994.00:.25:2019.75 ;
time06 = find(time == 2002) ;
time11 = find(time == 2012) ;
time05 = find(time == 2006) ;
time84 = find(time == 1984) ;

probsS(:,1) = (probsT(:,1) > 0.5)*1 + 0 ;
probsS(:,2) = (probsT(:,2) > 0.5)*1 + 0 ;

figure 
plot(time,(probsT(:,2)),'LineWidth',2)
axis tight
grid on
box on
% title("Estimated Probability of Bubbly Regime")
savefig('bubbly_regime')
exportgraphics(gcf,'graph/eps/bubbly_regime.eps')
exportgraphics(gcf,'graph/bubbly_regime.png','Resolution',300)

figure
subplot(2,1,1)
plot(time,At_mat(:,1),'LineWidth',2) ;
grid on
axis tight
title('生産性ショック')
subplot(2,1,2)
plot(time,-At_mat(:,2),'LineWidth',2) ;
title('選好ショック')
grid on
axis tight
% sgtitle("図表２: 各ショックの推計値")
exportgraphics(gcf,'graph/eps/shocks.eps')
exportgraphics(gcf,'graph/shocks.png','Resolution',300)

figure
plot(time,At_mat(:,1),'LineWidth',2) ;
grid on
axis tight
title('productivity shock')
savefig('productivity')
close gcf;

figure
plot(time,-At_mat(:,2),'LineWidth',2) ;
grid on
axis tight
title('preference shock')
savefig('preference')
close gcf;