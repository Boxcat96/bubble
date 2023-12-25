function [loglh, probsT, At_mat, At_mat2prob, bubble] = evallp_moments(phi,dds,gx,hx,nstates,ETAMATRIX,YY)

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

adjut = 1 ; % 1/10 ;
H       = 100*gx ;
ZZ1     = [H(1,:);H(31,:)] ;
ZZ2     = [H(2,:);H(32,:)*adjut] ;
ZC1     = H([1 5 7 9],:);
ZC2     = H([2 6 8 10],:) ;
RR      = 1e-7*eye(measerr); %zeros(measerr);          
RR(2,2) = 0.005 ;
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
Yt_mat  = zeros(nobs,4) ;

%% initialization
% random state draw
if rand < probs(1)
    ZZ  = ZZ1;
else
    ZZ  = ZZ2;
end
At = zeros(nstates,1);
Pt = zeros(nstates,nstates);
for tt=1:5
Pt = TT*Pt*TT' + RR*SIG*RR';
end

At = mvnrnd(At,Pt)';
% Condition on the first observation 
yhat = ZZ*At ;
nu   = YY(1,:) - yhat';
Ft   = ZZ*Pt*ZZ'+RR;
Ft   = 0.5*(Ft + Ft');
At   = At + (Pt*ZZ')/(Ft)*nu';
Pt   = Pt - (Pt*ZZ')/(Ft)*(Pt*ZZ')';

% compute likelihood with Kalman filter 

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

    for ii=1:2
        
    % Forecasting Step
    alphahat1  = TT*At1(ii,:)';
    alphahat2  = TT*At1(ii,:)';

    Phat1      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;
    Phat2      = TT*reshape(Pt1(ii,:),nstates,nstates)*TT' + SIG;

    yhat1      = DD1 + ZZ1*alphahat1 ;
    yhat2      = DD2 + ZZ2*alphahat2 ;

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
    At_mat2prob(t,:) = (1-probs(2))*(-0.1348) + 0.1348*probs(2)*(1+gx*At1(2,:)') ;
    At_mat(t,:)  = probs(1)*At1(1,:) + probs(2)*At1(2,:);
    Yt_mat(t,:)  = (ZC1*probs(1)*At1(1,:)' + ZC2*probs(2)*At1(2,:)')' ;
    probsT(t,:) = probs';
end

ssT_draw = [] ;

bubble = probsT(:,2).*(1+gx(end,:)*At_mat')' ;

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

