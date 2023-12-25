function [lik,spred] = linear_growth_smoother(gx,hx,nstates,~,~,ETAMATRIX,zz,tbreak,rhostruc)

% Objective: Compute smoothed shocks
% Growth model; linear approximation
% Pablo A. Guerron-Quintana
% Feb, 2 2013
% Revision Jan 2015
measerr = 2 ;
T   = size(zz,1) ;                                  % Number of time observations
R   = zeros(measerr) ;      % KF's measurement noise variance, no measure errors 
Q   = zeros(nstates,nstates) ;
Q(end-1:end,end-1:end) = ETAMATRIX(end-1:end,:) ;   % KF's process noise variance; see notes in pages 3
Q   = Q.^2 ;
S0  = zeros(nstates,1) ;                     % initial state
% Q   = [gx;eye(nstates)]*Q*[gx;eye(nstates)]' ;
% S0  = zeros(nstates+npsvar,1) ;                     % initial state

% format matrices

F       = hx ; % [zeros(size(gx,1),npsvar) gx*hx; zeros(size(hx,1),npsvar) hx] ;
H       = gx ;
H       = 100*H(1:8,:) ;

% Hfund   = [H(1,:);-H(3,:)+H(5,:)] ;
% Hbubb   = [H(2,:);-H(4,:)+H(6,:)] ;

Hfund   = [H(1,:);-H(7,:)+H(5,:)] ;
Hbubb   = [H(2,:);-H(8,:)+H(6,:)] ;

C       = zeros(nstates,1) ;
c       = zeros(measerr,1) ;  

% initialize lik

P0    = eye(nstates^2)/(eye(((nstates)^2))-kron(F,F))*reshape(Q,((nstates)^2),1);
P0    = reshape(P0,nstates,nstates);

lik     = 0 ;
sig     = P0 ;
shat    = S0 ;
ypred   = zeros(T,measerr) ;

for t=1:tbreak,    
  
  % Kalman Smoother step
  % =================================
  
  [shatnew,signew,sigpred(:,:,t),ypred(t,:)] = ksmoother(zz(t,:)',Hbubb,F,c,C,shat,sig,R,Q) ;

  shat = shatnew ;
  sig  = signew ;

  spred(t,:) = shatnew ;
  sigup(:,:,t) = signew ;
  
end;   % End of t loop.

for t=tbreak:T,    
  
  % Kalman Filter step
  % =================================
  
    [shatnew,signew,sigpred(:,:,t),ypred(t,:)] = ksmoother(zz(t,:)',Hfund,F,c,C,shat,sig,R,Q) ;

  shat = shatnew ;
  sig  = signew ;

  spred(t,:) = shatnew ;
  sigup(:,:,t) = signew ;
  
end;

ssmooth = zeros(size(spred)) ;
ssmooth(T,:) = spred(T,:) ;
Ftra = F' ;
sigsmooth = sigpred(:,:,T) ;

for t = T-1:-1:1
    isig = pinv(sigpred(:,:,t+1)) ; % (sigpred(:,:,t+1))\eye(size(sigpred(:,:,t+1))) ;
    ssmooth(t,:) = (spred(t,:)' + sigup(:,:,t)*Ftra*isig*(ssmooth(t+1,:)'-F*spred(t,:)'-C))' ;
    sigsmooth   = sigup(:,:,t) + sigup(:,:,t)*Ftra*isig*(sigsmooth - sigpred(:,:,t+1))*(isig)'*F*sigup(:,:,t)' ;   
end

% Recover shocks
rhos = [rhostruc.rhoa rhostruc.rhob] ;

% order in shocks is reversed
% b a
shocks = zeros(measerr,T) ;

for t = 1:T
    ysmooth(t,:) = gx*ssmooth(t,:)' ;
    cybsmooth(t,:) = Hbubb*ssmooth(t,:)' ;
    cyfsmooth(t,:) = Hfund*ssmooth(t,:)' ;
    for ix = 1:measerr
        if t == 1
            shocks(ix,t) = ssmooth(t,end-ix+1) ;
        else
            shocks(ix,t) = ssmooth(t,end-ix+1) - rhos(end-ix+1)*ssmooth(t-1,end-ix+1) ;
        end
    end
end

time = 1947.25:.25:2016.75 ;
figure(33)
% hold on
%plot(time,[cybsmooth(1:tbreak,1); cyfsmooth(tbreak+1:end,1)],'LineWidth',2)
%plot(time,cybsmooth(:,1),'--','LineWidth',2)
% hold off
% plot(time,[cybsmooth(1:tbreak,2); cyfsmooth(tbreak+1:end,2)],'LineWidth',2)

% change mbt value!!!! depends on each estimation
mbts = 100*0.0204 ; % 100 x mbs
[mbtrend,mbcycle] = hpfilter(mbts*(1+ysmooth(1:tbreak,end)/100)) ;
hold on
plot(time(1:tbreak),mbts*(1+ysmooth(1:tbreak,end)/100),'LineWidth',2)
plot(time(1:tbreak),mbtrend,'r--','LineWidth',2)
legend('Real Bubble','HP-Trend')
xlim([1947,1985])
ylim([1.7,2.3])
hold off
 
fprintf('Relative Std real bubble %f ',100*std(mbcycle)/std(zz(:,1))) % x 100 to make it compatible with output growth units

[std(shocks(1,1:tbreak)) std(shocks(2,1:tbreak))]
[std(shocks(1,tbreak+1:end)) std(shocks(2,tbreak+1:end))]
[std(shocks(1,1:tbreak))/std(shocks(1,tbreak+1:end)) std(shocks(2,1:tbreak))/std(shocks(2,tbreak+1:end))]

figure(44)
subplot(2,2,1)
plot(time,zz(:,1),'LineWidth',2)
title('Output Growth')
axis tight
subplot(2,2,2)
plot(time,zz(:,2),'LineWidth',2)
title('Consumption-Investment Ratio')
axis tight
subplot(2,2,3)
title('Shocks Pre-1984')
hold on
yyaxis left
plot(time(1:tbreak), shocks(1,1:tbreak),'LineWidth',2)
ylim([-6,6])
yyaxis right
plot(time(1:tbreak), shocks(2,1:tbreak),'--','LineWidth',2)
ylim([-1.5,1.5])
% axis tight
legend('Supply Shock -- Left Axis','Demand Shock -- Right Axis')
hold off
xlim([1947,1985])
subplot(2,2,4)
title('Shocks Post-1984')
hold on 
yyaxis left
plot(time(tbreak+1:end), shocks(1,tbreak+1:end),'LineWidth',2)
ylim([-4,4])
yyaxis right
plot(time(tbreak+1:end), shocks(2,tbreak+1:end),'--','LineWidth',2)
ylim([-1,1])
xlim([1984,2017])
% axis tight
hold off

allshocks = 1 ;

if allshocks == 1
    
    load('../estimation_bubble/shocks_bubble.mat') ;
    load('../estimation_nobubble/shocks_nobubble.mat')
    
    figure(55)
    subplot(2,1,2)
    hold on
    plot(time,shocks(1,:),'LineWidth',2)
    plot(time,shocks_bubble(1,:),'--','LineWidth',2)
    plot(time,shocks_nobubble(1,:),':','LineWidth',2)
    hold off
    axis tight
    title('Productivity Shocks')
    subplot(2,1,1)
    hold on
    plot(time,shocks(2,:),'LineWidth',2)
    plot(time,shocks_bubble(2,:),'--','LineWidth',2)
    plot(time,shocks_nobubble(2,:),':','LineWidth',2)
    hold off
    axis tight
    legend('Benchmark','Permanent Bubble','Bubbleless')
    title('Preference Shocks')
end    
    