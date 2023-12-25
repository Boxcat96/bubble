function [lik,ypred] = linear_growth(gx,hx,nstates,~,~,ETAMATRIX,tbreak,zz)

% Objective: Compute likelihood function
% Growth model; linear approximation
% Pablo A. Guerron-Quintana
% August, 25 2013
% Revision Dec 2014
measerr = 2 ; 
T   = size(zz,1) ;                                  % Number of time observations
R   = zeros(measerr) ; % zeros(measerr) ;      % KF's measurement noise variance, no measure errors 
Q   = zeros(nstates,nstates) ;
Q(end-1:end,end-1:end) = ETAMATRIX(end-1:end,:) ;   % KF's process noise variance; see notes in pages 3
Q   = Q.^2 ;
S0  = zeros(nstates,1) ;

% format matrices

F       = hx ;
H       = gx ;
H       = 100*H(1:8,:) ;

% Hfund   = [H(1,:);-H(3,:)+H(5,:)] ;
% Hbubb   = [H(2,:);-H(4,:)+H(6,:)] ;

Hfund   = [H(1,:);-H(7,:)+H(5,:)] ;
Hbubb   = [H(2,:);-H(8,:)+H(6,:)] ;

C       = zeros(nstates,1) ; % no constant in state equation
c       = zeros(measerr,1) ;        % 

% initialize lik

% P0    = inv(eye(((nstates)^2))-kron(F,F))*reshape(Q,((nstates)^2),1);
P0    = (eye(nstates^2)/(eye(((nstates)^2))-kron(F,F)))*reshape(Q,((nstates)^2),1);
P0    = reshape(P0,nstates,nstates);
% save P0 P0
% load P0

lik     = 0 ;
sig     = P0 ;
shat    = S0 ;
ypred   = zeros(T,measerr) ;

for t=1:tbreak    
  
  % Kalman Filter step
  % =================================
  
  [shatnew,signew,liko,ypred(t,:)]=kfilteralt(zz(t,:)',Hbubb,F,c,C,shat,sig,R,Q,measerr) ;

  shat = shatnew ;
  sig  = signew ;
  lik = lik + liko ;
  
end;   % End of t loop.

for t=tbreak:T    
  
  % Kalman Filter step
  % =================================
  
  [shatnew,signew,liko,ypred(t,:)]=kfilteralt(zz(t,:)',Hfund,F,c,C,shat,sig,R,Q,measerr) ;

  shat = shatnew ;
  sig  = signew ;
  lik = lik + liko ;
  
end;
