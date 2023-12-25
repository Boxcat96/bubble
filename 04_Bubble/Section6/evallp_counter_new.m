function [loglh, probsT, At_mat, At_mat2prob, bubble] = evallp_counter_new(phi,dds,gx,hx,nstates,ETAMATRIX,YY)

At_mat2prob = [] ;
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
DD1(1,1) = dds(1) ; 
DD1(2,1) = dds(2) ; 
DD2     = zeros(nstates,1) ;
DD2(1,1) = dds(3) ; 
DD2(2,1) = dds(4) ; 

H       = 100*gx ;
ZZ1     = [H(3,:);H(31,:)] ;
ZZ2     = [H(4,:);H(32,:)] ; 
ZZy1    = [H(3,:);0 0] ;
ZZy2    = [H(4,:);0 0] ;
ZZg1    = [H(1,:);0 0] ;
ZZg2    = [H(2,:);0 0] ;
RR      = 1e-7*eye(measerr);          
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
ynoshk1c = ynoshk ; 
load steadystate ;
ybt = 100*(ybt) ; yft = 100*(yft) ;
gstlag = DD1(1) ;
yslag  = yft ;

probby = 0.6 ;

for t=1:nobs
    
    yhat1c(:,t)      = DD1 + ZZ1*At_mat(t,:)' - ylag(:,1) + ghat(:,1) ;
    yhat2c(:,t)      = DD2 + ZZ2*At_mat(t,:)' - ylag(:,2) + ghat(:,2) ;
    ylag(:,1)        = ZZy1*At_mat(t,:)';
    ylag(:,2)        = ZZy2*At_mat(t,:)';
    ghat(:,1)        = ZZg1*At_mat(t,:)' ;
    ghat(:,2)        = ZZg2*At_mat(t,:)' ;
    
    ynoshk(t)        = yft*((1 - probsT(t,2))>probby) + ybt*(probsT(t,2)>1-probby) + gstlag - yslag ;
    gstlag           = DD1(1)*((1 - probsT(t,2))>probby) + DD2(1)*(probsT(t,2)>1-probby) ;
    yslag            = yft*((1 - probsT(t,2))>probby) + ybt*(probsT(t,2)>1-probby) ;
    ynoshk1c(t)      = DD1(1) ;
    
end

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

% Growth in steady state fundamental regime in bubbly equilibrium
gfbss = gft ;

% Load steady state fundatamental model
load steadystatefund ;
% Growth in steady state in fundamental model
gfss = gft ;
ybt = 100*(ybt) ; yft = 100*(yft) ;
gstlag = DD1(1) ;
yslag  = yft ;

for t=1:nobs
    
    yhat1cp(:,t)      = DD1 + 100*(gfss - gfbss) ; % DD1 + ZZ1*At_mat(t,:)' ;
    
end

%% data_start
time1 = 1994.00:.25:2019.75 ;
time97 = find(time1 == 1996) ;
time05 = find(time1 == 2004) ;
time10 = find(time1 == 2012) ;

counterf97 = zeros(3,length(time1(time97:time05))) ;
counterf97(:,1) = 1 ;
counterf = zeros(3,length(time1(time05:time10))) ;
counterf(:,1) = 1 ;

for tt = 1:length(time1(time97:time05))
    counterf97(1,tt+1) = counterf97(1,tt)*(1+yhat1cp(1,time97+tt-1)/100) ;
    counterf97(2,tt+1) = counterf97(2,tt)*(1+ynoshk1c(1,time97+tt-1)/100) ;
    counterf97(3,tt+1) = counterf97(3,tt)*(1+ynoshk(1,time97+tt-1)/100) ;
end

for tt = 1:length(time1(time05:time10))
    counterf(1,tt+1) = counterf(1,tt)*(1+yhat1cp(1,time05+tt-1)/100) ;
    counterf(2,tt+1) = counterf(2,tt)*(1+ynoshk1c(1,time05+tt-1)/100) ;
    counterf(3,tt+1) = counterf(3,tt)*(1+ynoshk(1,time05+tt-1)/100) ;
end

% subplot(3,1,3)

% figure('Units','centimeters','Position',[1 1 40 40])
figure
sgtitle("Figure 10: Counterfactual Simulations")
subplot(2,1,2)
hold on
plot(time1(time05:time10),100*log(counterf(1,1:end-1)),'k:','LineWidth',3)
plot(time1(time05:time10),100*log(counterf(2,1:end-1)),'r--','LineWidth',3)
plot(time1(time05:time10),100*log(counterf(3,1:end-1)),'b','LineWidth',3)
hold off
grid on
axis tight
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','northwest','FontSize',13)
title('Housing bubble')

subplot(2,1,1)
hold on
plot(time1(time97:time05),100*log(counterf97(1,1:end-1)),'k:','LineWidth',3)
plot(time1(time97:time05),100*log(counterf97(2,1:end-1)),'r--','LineWidth',3)
plot(time1(time97:time05),100*log(counterf97(3,1:end-1)),'b','LineWidth',3)
hold off
grid on
axis tight
set(gca,'FontSize',12)
box on
legend({'fundamental','recurrent-no bubble','benchmark'},'Location','northwest','FontSize',13)
title('"After-Bubble" bubble')

exportgraphics(gcf,'graph/eps/counterfactual.eps')
exportgraphics(gcf,'graph/counterfactual.png','Resolution',300)

return

