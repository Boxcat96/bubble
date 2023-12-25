function [loglh, probsT, At_mat, At_mat2prob, bubble] = evallp_counter(phi,dds,gx,hx,nstates,ETAMATRIX,YY)

close all

%% house cleaning
% YY is data: zz in your notation
% phi is a 2-dimensional vector: markov transition probabilities

%% state-space representation
% Y_t = DD + ZZ(S_t)*a_t + u_t, u_t~N(0,RR) 
% a_t = TT*a_{t-1} + e_t,      e_t~N(0,SIG)
% S_t = {1,2}

% This version allows for break in volatilities in 1984

[nobs,measerr] = size(YY);   
% measurement equation  
DD1     = zeros(nstates,1) ;
DD1(1,1) = dds(1) ; % 0.01 ;
DD1(2,1) = dds(2) ; %-.025 ;
DD2     = zeros(nstates,1) ;
DD2(1,1) = dds(3) ; % 0.5 ;
DD2(2,1) = dds(4) ; % 0.05 ;

H       = 100*gx ;
ZZ1     = [H(3,:);H(31,:)] ;
ZZ2     = [H(4,:);H(32,:)] ; 
ZZy1    = [H(3,:);0 0] ;
ZZy2    = [H(4,:);0 0] ;
ZZg1    = [H(1,:);0 0] ;
ZZg2    = [H(2,:);0 0] ;
RR      = 1e-7*eye(measerr); %zeros(measerr);          
RR(2,2) = 1 ; % 5\% meas error
% transition equation
TT      = hx ;
SIG     = zeros(nstates,nstates);
SIG(end-1:end,end-1:end) = ETAMATRIX(end-1:end,:);   
SIG     = SIG.^2;

%% Markov Transition Matrix
strans = [phi(1) 1-phi(2) ; 1-phi(1) phi(2)];
probs  = [(1-phi(2))/(2-phi(1)-phi(2));(1-phi(1))/(2-phi(1)-phi(2))];

% likelihood at t=0
loglh  = 0; 

% store
At_mat1 = zeros(nobs,nstates);
At_mat2 = zeros(nobs,nstates);
%At_mat2prob = zeros(size(hx)) ;
At_mat  = zeros(nobs,nstates);

%% initialization
% random state draw
if rand < probs(1)
    ZZ  = ZZ1;
    ZZy = ZZy1 ;
    ZZg = ZZg1 ;
    DD  = DD1 ;
else
    ZZ  = ZZ2;
    ZZy = ZZy2 ;
    ZZg = ZZg2;
    DD  = DD2 ;
end
At = zeros(nstates,1);
Pt = zeros(nstates,nstates);
for tt=1:5
    Pt = TT*Pt*TT' + RR*SIG*RR';
end

At = mvnrnd(At,Pt)';
% Condition on the first observation
ghatt = ZZg*At ;
ylagt= ZZy*At ;
yhat = DD + ZZ*At - ylagt + ghatt ;
nu   = YY(1,:) - yhat';
Ft   = ZZ*Pt*ZZ'+RR;
Ft   = 0.5*(Ft + Ft');
At   = At + (Pt*ZZ')/(Ft)*nu';
Pt   = Pt - (Pt*ZZ')/(Ft)*(Pt*ZZ')';

% compute likelihood with Kalman filter 

ylag(:,1) = ylagt ;
ylag(:,2) = ylagt ;
ghat(:,1) = ghatt ;
ghat(:,2) = ghatt ;
At1 = repmat(At',2,1);
Pt1 = repmat(Pt(:)',2,1);

% Initialize the matrix that store the info from the forward simulation   
probsT      = zeros(nobs,2);
probsT(1,:) = probs';     %The first element is unconditional prob

for t=1:nobs

    % Forecasting the state probabilities
    probx  = kron(probs,ones(2,1));
    probx  = strans(:).*probx;

    loglht   = zeros(4,1);
    At       = zeros(4,nstates);
    Pt       = zeros(4,nstates^2);
    
    % Loop over t-1 regimes
    for ii=1:2
        
        % Forecasting Step
        alphahat1  = TT*At1(ii,:)';
        alphahat2  = TT*At1(ii,:)';

        Phat1      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;
        Phat2      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;

        % ghat = ZZg*At ;
        % ylag = ZZy*At ;
        % yhat = ZZ*At - ylag + ghat ;
        yhat1      = DD1 + ZZ1*alphahat1 - ylag(:,ii) + ghat(:,ii) ;
        yhat2      = DD2 + ZZ2*alphahat2 - ylag(:,ii) + ghat(:,ii) ;

        nu1        = YY(t,:) - yhat1';
        nu2        = YY(t,:) - yhat2';

        Ft1        = ZZ1*Phat1*ZZ1' + RR;
        Ft2        = ZZ2*Phat2*ZZ2' + RR;

        Ft1        = 0.5*(Ft1 + Ft1');
        Ft2        = 0.5*(Ft2 + Ft2');


        % Compute likelihood
        iFt1       = eye(size(Ft1))/(Ft1) ;
        iFt2       = eye(size(Ft2))/(Ft2) ;
        loglht(2*(ii-1)+1,1) = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft1))) - 0.5*nu1*iFt1*nu1' );
        loglht(2*ii,1)       = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft2))) - 0.5*nu2*iFt2*nu2' );

        % Updating Step
        At(2*(ii-1)+1,:)   = ( alphahat1 + Phat1*ZZ1'*(iFt1)*nu1' )';
        At(2*ii,:)         = ( alphahat2 + Phat2*ZZ2'*(iFt2)*nu2' )';

        Pt(2*(ii-1)+1,:)   = vec(Phat1 - Phat1*ZZ1'*(iFt1)*ZZ1*Phat1' )';
        Pt(2*ii,:)         = vec(Phat2 - Phat2*ZZ2'*(iFt2)*ZZ2*Phat2' )';

    end

    % loglh computation
    loglhtmax = max(loglht);
    loglh     = loglh + loglhtmax + log(probx'*exp(loglht-loglhtmax));

    % Updating mixture probabilities and the state probabilities
    probx = (probx.*exp(loglht-loglhtmax)) / ( probx' * exp(loglht-loglhtmax) );
    probs = [probx(1) + probx(3); probx(2) + probx(4)];   

    for ii=1:2
    At1(ii,:) = (probx(ii)/probs(ii)) * At(ii,:)+  (probx(ii+2)/probs(ii)) * At(ii+2,:);
    Pt1(ii,:) = vec( (probx(ii)/probs(ii))*(reshape(Pt(ii,:),nstates,nstates)+At(ii,:)'*At(ii,:))...
              + (probx(ii+2)/probs(ii))*(reshape(Pt(ii+2,:),nstates,nstates)+At(ii+2,:)'*At(ii+2,:))...
              - (At1(ii,:)'*At1(ii,:))');
    end      

    At_mat1(t,:) = At1(1,:);
    At_mat2(t,:) = At1(2,:);
    ylag(:,1)    = ZZy1*At1(1,:)' ;
    ylag(:,2)    = ZZy2*At1(2,:)' ;
    ghat(:,1)    = ZZg1*At1(1,:)' ;
    ghat(:,2)    = ZZg2*At1(2,:)' ;
    At_mat2prob(t,:) = (1-probs(2))*(-0.1348) + 0.1348*probs(2)*(1+gx*At1(2,:)') ;
    At_mat(t,:)  = probs(1)*At1(1,:) + probs(2)*At1(2,:);

    probsT(t,:) = probs';
end

bubble = probsT(:,2).*(1+gx(end,:)*At_mat')' ;

% Now, counterfactuals

% 1. All time in one regime but with expectation of switching
yhat1c = zeros(2,nobs) ;
yhat2c = yhat1c ;
ylag(:,1) = ylagt ;
ylag(:,2) = ylagt ;
ghat(:,1) = ghatt ;
ghat(:,2) = ghatt ;

ynoshk  = zeros(1,nobs) ;

load steadystate ;
ybt = 100*(ybt) ; yft = 100*(yft) ;
gstlag = DD1(1) ;
yslag  = yft ;

for t=1:nobs
    
    yhat1c(:,t)      = DD1 + ZZ1*At_mat(t,:)' - ylag(:,1) + ghat(:,1) ;
    yhat2c(:,t)      = DD2 + ZZ2*At_mat(t,:)' - ylag(:,2) + ghat(:,2) ;
    ylag(:,1)        = ZZy1*At_mat(t,:)';
    ylag(:,2)        = ZZy2*At_mat(t,:)';
    ghat(:,1)        = ZZg1*At_mat(t,:)' ;
    ghat(:,2)        = ZZg2*At_mat(t,:)' ;
    
    ynoshk(t)        = yft*((1 - probsT(t,2))>0.5) + ybt*(probsT(t,2)>0.5) + gstlag - yslag ;
    gstlag           = DD1(1)*((1 - probsT(t,2))>0.5) + DD2(1)*(probsT(t,2)>0.5) ;
    yslag            = yft*((1 - probsT(t,2))>0.5) + ybt*(probsT(t,2)>0.5) ;
    
end

time1 = 1984:.25:2017.75 ;

figure
subplot(2,1,1)
plot(time1,ynoshk,'LineWidth',2) ;
axis tight
subplot(2,1,2)
plot(time1,cumsum(ynoshk),'LineWidth',2) ;
axis tight
pause

baseline = yhat1c(1,:).*(1 - probsT(:,2)') + yhat2c(1,:).*probsT(:,2)' ;

% hold on
% plot(time1,yhat1c(1,:),'--','LineWidth',2)
% plot(time1,yhat2c(1,:),'r','LineWidth',2)
% plot(time1,YY(:,1),'k','LineWidth',2)
% legend('fundamental','bubbly','data')
% hold off

[mean(yhat1c(1,:)) mean(yhat2c(1,:)) mean(YY(:,1)) mean(baseline)]
[std(yhat1c(1,:)) std(yhat2c(1,:)) std(YY(:,1)) std(baseline)]

% 2. Permanent regime

load gxhxcounter
H       = 100*gxc ;
ZZ1     = [H(3,:);H(31,:)] ;
ZZ2     = [H(4,:);H(32,:)] ; 
ZZy1    = [H(3,:);0 0] ;
ZZy2    = [H(4,:);0 0] ;
ZZg1    = [H(1,:);0 0] ;
ZZg2    = [H(2,:);0 0] ;

yhat1cp  = zeros(2,nobs) ;
yhat2cp  = yhat1cp ;
ylagc    = ylagt ;
ghatc    = ghatt ;

% from compute_steady_state_growth file
% gf = 0.0056 at fundamental but bubbly eq
% gf = 0.0109 at fundamental

for t=1:nobs
    
    yhat1cp(:,t)      = DD1 + 100*(0.0109 - 0.0056) + ZZ1*At_mat(t,:)' - ylagc + ghatc ; % DD1 + ZZ1*At_mat(t,:)' ;
    yhat2cp(:,t)      = DD2 + ZZ2*At_mat(t,:)' ;
    ylagc             = ZZy1*At_mat(t,:)' ;
    ghatc             = ZZg1*At_mat(t,:)' ;
    
end

time1 = 1984:.25:2017.75 ;
time05 = find(time1 == 2005) ;
time10 = find(time1 == 2011.75) ;
time10 = find(time1 == 2017.5) ;
% hold on
% plot(time1,yhat1cp(1,:),'--','LineWidth',2)
% plot(time1,yhat2cp(1,:),'r','LineWidth',2)
% plot(time1,YY(:,1),'k','LineWidth',2)
% legend('fundamental-per','bubbly-per')
% hold off

[mean(yhat1cp(1,:)) mean(yhat2cp(1,:)) mean(YY(:,1))]
[std(yhat1cp(1,:)) std(yhat2cp(1,:)) std(YY(:,1))]

figure('Units','centimeters','Position',[1 1 30 30])

subplot(3,1,1)
hold on
% plot(time1,yhat1cp(1,:)-yhat1c(1,:),'--','LineWidth',2)
plot(time1,yhat1cp(1,:),'k:','LineWidth',3) ;
plot(time1,yhat1c(1,:),'r--','LineWidth',3)
plot(time1,baseline,'LineWidth',3)
hold off
axis tight
grid on
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','southwest','FontSize',13)
%p = get(gca, 'Position')
%h = axes('Parent', gcf, 'Position', [p(1)+.2 p(2)+.05 p(3)-.5 p(4)-.5]);
subplot(3,1,2)
hold on
plot(time1(time05:time10),yhat1cp(1,time05:time10),'k:','LineWidth',3)
plot(time1(time05:time10),yhat1c(1,time05:time10),'r--','LineWidth',3)
plot(time1(time05:time10),baseline(1,time05:time10),'LineWidth',3)
hold off
grid on
axis tight
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','southwest','FontSize',13)

% figure 
% 
% hold on
% plot(time1(time05:time10),cumsum(yhat1cp(1,time05:time10)),'k:','LineWidth',3)
% plot(time1(time05:time10),cumsum(yhat1c(1,time05:time10)),'r--','LineWidth',3)
% plot(time1(time05:time10),cumsum(baseline(1,time05:time10)),'LineWidth',3)
% hold off
% grid on
% axis tight
% legend({'fundamental-per','fundament-switch','baseline'},'Location','southwest','FontSize',13)


counterf = zeros(3,length(time1(time05:time10))) ;
counterf(:,1) = 1 ;

for tt = 1:length(time1(time05:time10))
    counterf(1,tt+1) = counterf(1,tt)*(1+yhat1cp(1,time05+tt-1)/100) ;
    counterf(2,tt+1) = counterf(2,tt)*(1+yhat1c(1,time05+tt-1)/100) ;
    counterf(3,tt+1) = counterf(3,tt)*(1+baseline(1,time05+tt-1)/100) ;
end

subplot(3,1,3)

hold on
plot(time1(time05:time10),counterf(1,1:end-1),'k:','LineWidth',3)
plot(time1(time05:time10),counterf(2,1:end-1),'r--','LineWidth',3)
plot(time1(time05:time10),counterf(3,1:end-1),'LineWidth',3)
hold off
grid on
axis tight
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','southwest','FontSize',13)

figure('Units','centimeters','Position',[1 1 30 30])

hold on
% plot(time1,yhat1cp(1,:)-yhat1c(1,:),'--','LineWidth',2)
plot(time1,yhat1cp(1,:),'k:','LineWidth',3) ;
plot(time1,yhat1c(1,:),'r-x','LineWidth',3)
plot(time1,baseline,'LineWidth',3)
hold off
axis tight
grid on
set(gca,'FontSize',12)
box on
legend({'no chance of bubbles','no bubbles by chance','benchmark'},'Location','southwest','FontSize',13)
%p = get(gca, 'Position')
%h = axes('Parent', gcf, 'Position', [p(1)+.2 p(2)+.05 p(3)-.5 p(4)-.5]);
title('GDP Growth')

pause 

figure('Units','centimeters','Position',[1 1 30 50])

subplot(2,1,1)
hold on
% plot(time1,yhat1cp(1,:)-yhat1c(1,:),'--','LineWidth',2)
plot(time1,yhat1cp(1,:),'k:','LineWidth',3) ;
plot(time1,yhat1c(1,:),'r--','LineWidth',3)
plot(time1,baseline,'LineWidth',3)
hold off
axis tight
grid on
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','southwest','FontSize',13)
%p = get(gca, 'Position')
%h = axes('Parent', gcf, 'Position', [p(1)+.2 p(2)+.05 p(3)-.5 p(4)-.5]);
title('GDP Growth')

subplot(2,1,2)

hold on
plot(time1(time05:time10),counterf(1,1:end-1),'k:','LineWidth',3)
plot(time1(time05:time10),counterf(2,1:end-1),'r--','LineWidth',3)
plot(time1(time05:time10),counterf(3,1:end-1),'LineWidth',3)
hold off
grid on
axis tight
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','northwest','FontSize',13)
title('GDP Level') 

return

figure

hold on
plot(time1,yhat2cp(1,:) - yhat2c(1,:),'LineWidth',2)
% plot(time1,yhat2c(1,:),'r','LineWidth',2)
% plot(time1,YY(:,1),'k','LineWidth',2)
% legend('bubbly-per','bubbly-switch')
hold off
axis tight; grid on
legend('Bubbly Permament vs Blind')

figure

hold on
plot(time1,yhat1cp(1,:) - yhat1c(1,:),'LineWidth',2)
hold off
axis tight; grid on
return

% size(At_mat2prob)

% for t=148:nobs
% 
%     % Forecasting the state probabilities
%     probx  = kron(probs,ones(2,1));
%     probx  = strans(:).*probx;
% 
%     loglht   = zeros(4,1);
%     At       = zeros(4,nstates);
%     Pt       = zeros(4,nstates^2);
% 
%     for ii=1:2
%         
%     % Forecasting Step
%     alphahat1  = TT*At1(ii,:)';
%     alphahat2  = TT*At1(ii,:)';
% 
%     Phat1      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;
%     Phat2      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;
% 
%     yhat1      = DD1 + ZZ1*alphahat1 ;
%     yhat2      = DD2 + ZZ2*alphahat2 ;
% 
%     nu1        = YY(t,:) - yhat1';
%     nu2        = YY(t,:) - yhat2';
% 
%     Ft1        = ZZ1*Phat1*ZZ1' + RR;
%     Ft2        = ZZ2*Phat2*ZZ2' + RR;
% 
%     Ft1        = 0.5*(Ft1 + Ft1');
%     Ft2        = 0.5*(Ft2 + Ft2');
% 
% 
%     % Compute likelihood
%     iFt1       = eye(size(Ft1))/(Ft1) ;
%     iFt2       = eye(size(Ft2))/(Ft2) ;
%     loglht(2*(ii-1)+1,1) = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft1))) - 0.5*nu1*iFt1*nu1' );
%     loglht(2*ii,1)       = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft2))) - 0.5*nu2*iFt2*nu2' );
% 
%     % Updating Step
%     At(2*(ii-1)+1,:)   = ( alphahat1 + Phat1*ZZ1'*(iFt1)*nu1' )';
%     At(2*ii,:)         = ( alphahat2 + Phat2*ZZ2'*(iFt2)*nu2' )';
% 
%     Pt(2*(ii-1)+1,:)   = vec(Phat1 - Phat1*ZZ1'*(iFt1)*ZZ1*Phat1' )';
%     Pt(2*ii,:)         = vec(Phat2 - Phat2*ZZ2'*(iFt2)*ZZ2*Phat2' )';
% 
%     end
% 
%     % loglh computation
%     loglhtmax = max(loglht);
%     loglh     = loglh + loglhtmax + log(probx'*exp(loglht-loglhtmax));
% 
%     % Updating mixture probabilities and the state probabilities
%     probx = (probx.*exp(loglht-loglhtmax)) / ( probx' * exp(loglht-loglhtmax) );
%     probs = [probx(1) + probx(3); probx(2) + probx(4)];   
% 
%     for ii=1:2
%     At1(ii,:) = (probx(ii)/probs(ii)) * At(ii,:)+  (probx(ii+2)/probs(ii)) * At(ii+2,:);
%     Pt1(ii,:) = vec( (probx(ii)/probs(ii))*(reshape(Pt(ii,:),nstates,nstates)+At(ii,:)'*At(ii,:))...
%               + (probx(ii+2)/probs(ii))*(reshape(Pt(ii+2,:),nstates,nstates)+At(ii+2,:)'*At(ii+2,:))...
%               - (At1(ii,:)'*At1(ii,:))');
%     end      
% 
%     At_mat1(t,:) = At1(1,:);
%     At_mat2(t,:) = At1(2,:);
%     At_mat(t,:)  = probs(1)*At1(1,:) + probs(2)*At1(2,:);
% 
%     probsT(t,:) = probs';
% end

ssT_draw = [] ;

% initialize rand
% potential seeds
% seeds = [10 12 74 85] ;
seeds = 10 ;


for ix = 1:length(seeds)

    rng(seeds(ix))

% Smoother
 t = nobs; 
   probs1     = probsT(t,1);
   draws1     = rand(1,1) <= probs1;
   ssT_draw(t,:) = [draws1  (1-draws1)];
   
   while t ~= 1;
  
      t = t-1;
      probs1 = ssT_draw(t+1,:)*strans(:,1)*probsT(t,1)/( ssT_draw(t+1,:)*strans*probsT(t,:)' );
      draws1 = rand(1,1) <= probs1;
      ssT_draw(t,:) = [draws1  (1-draws1)];
   
   
   end   

%  figure(ix)  
%    
%  hold on
%  time = 1947.25:.25:2016.75 ;
%  plot(time,probsT(:,2),'LineWidth',2)
%  plot(time,ssT_draw(:,2),'LineWidth',2)
%  axis tight
%  hold off
%  
%  pause
 % close
 
end
  
end

