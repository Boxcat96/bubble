% This procedure simulates from the latent vector of regimes and computes 
% the likelihood of the draws
% loglh: log likelihood
% retcode: determinacy or indeterminacy

function [ loglh, probsT, At_mat, At_mat2prob, bubble ] = evallp_mod1(phi,dds,gx,hx,nstates,ETAMATRIX,YY)

[nobs,measerr] = size(YY);   
% measurement equation
% 1 & 2 bubble, fundamental
% DD1     = zeros(nstates,1) ;
% DD1(1,1) = 0.5 ;
% DD1(2,1) = -0.1 ;
% DD2     = zeros(nstates,1) ;
% DD2(1,1) = 0.9 ;
% DD2(2,1) = 0.6 ;
% % 3 pure fundamental
% DD3     = zeros(nstates,1) ;
% DD3(1,1) = 1.0 ;
% DD3(2,1) = -0.1 ;

% DD1     = zeros(nstates,1) ;
% DD1(1,1) = 0.25 ;
% DD1(2,1) = 1.0 ;
% DD2     = zeros(nstates,1) ;
% DD2(1,1) = 0.9 ;
% DD2(2,1) = 1.3 ;
% % 3 pure fundamental
% DD3     = zeros(nstates,1) ;
% DD3(1,1) = 1.0 ;
% DD3(2,1) = 0.5 ;

DD1     = zeros(nstates,1) ;
DD1(1,1) = dds(1) ; %0.35 ;
DD1(2,1) = dds(2) ; %-0.05 ;
DD2     = zeros(nstates,1) ;
DD2(1,1) = dds(3) ; %0.5 ;
DD2(2,1) = dds(4) ; % 0.15 ;
% 3 pure fundamental
DD3     = zeros(nstates,1) ;
DD3(1,1) = dds(5) ; %0.8 ;
DD3(2,1) = dds(6) ; %0.00 ;

adjut1  = 1 ; %1/100 ;
adjut2  = 1 ;
H       = 100*gx ;
ZZ1     = [H(1,:);H(31,:)] ;
ZZ2     = [H(2,:);H(32,:)*adjut1] ;
ZZ3     = adjut2*[H(1,:);H(31,:)] ;
RR      = 1e-8*eye(measerr); %zeros(measerr);  
RR(2,2) = 0.001 ;

% transition equation
TT      = hx ;
SIG     = zeros(nstates,nstates);
SIG(end-1:end,end-1:end) = ETAMATRIX(end-1:end,:);   
SIG     = SIG.^2;

loglh     = 0;
loglhzero = -1E8;
probsT    = zeros(nobs,3);
 
%Markov Transition Matrix
prob3s = 1/2 ;
strans  = [prob3s*phi(1) prob3s*(1-phi(2)) prob3s/2; ...
          prob3s*(1-phi(1)) prob3s*phi(2) prob3s/2; ...
          1-prob3s      1-prob3s          1-prob3s] ;
          
probs = ss_prob(strans) ;
           
%% Solve the DSGE model 
%  Y_t = ZZ*X_t
%  X_t = TT*X_t-1+RR*e_t   e_t~N(0,SIG)  

% KK   = 0*eye(nobs); KK(end,end) = .0005;
% HH = KK;

r_rand = rand;
index = find((r_rand-cumsum(probs)).^2==min(r_rand-cumsum(probs)).^2);
ZZ  = ZZ1*(index == 1) + ZZ2*(index == 2) + ZZ3*(index == 3) ;
DD  = DD1*(index == 1) + DD2*(index == 2) + DD3*(index == 3) ;
        
nstate = size(TT,1);
At = zeros(nstate,1);
At_mat  = zeros(nobs,nstates);
%Pt = reshape(inv(eye(nstate^2)-kron(TT,TT))*reshape(RR*SIG*RR',nstate^2,1),nstate,nstate);
Pt = zeros(nstate,nstate);
for ii=1:5
    Pt = TT*Pt*TT' + SIG;
end

% Condition on the first observation 
yhat = DD + ZZ*At ;
nu   = YY(1,:) - yhat';
Ft   = ZZ*Pt*ZZ' + RR ;
Ft   = 0.5*(Ft + Ft');
At   = At + (Pt*ZZ')*inv(Ft)*nu';
Pt   = Pt - (Pt*ZZ')*inv(Ft)*(Pt*ZZ')';    

% compute likelihood with Kalman filter 

At1 = At';
Pt1 = vec(Pt)';
probx = 1;
     
% Initialize the matrix that store the info from the forward simulation
   
probsT = zeros(nobs,3);
probsT(1,:) = probs';     %The first element is unconditional prob
 
for t=1:nobs

%        if t==here
%            ZZ = ZZ_;
%            YY = YY_;
%            HH = KK;
%        end
      % Forecasting the state probabilities
      
    nx = size(At1,1);
    probs = strans*probs;
    probx2 = kron(probs,probx);

    loglht   = zeros(2*nx,1);
    At       = zeros(2*nx,nstate);
    Pt       = zeros(2*nx,nstate^2);
    yhat     = zeros(size(YY,2),1);
    Ft       = zeros(size(YY,2),size(YY,2));
      
        for i=1:nx
                    
            % Forecasting Step
               
            alphahat1  = TT*At1(i,:)';
            % alphahat2  = TT2*At1(i,:)';
            % alphahat3  = TT3*At1(i,:)';
            
            Phat1      = TT*reshape(Pt1(i,:),nstate,nstate)*TT' + SIG;
            %Phat2      = TT2*reshape(Pt1(i,:),nstate,nstate)*TT2' + SIG;
            %Phat3      = TT3*reshape(Pt1(i,:),nstate,nstate)*TT3' + SIG;
            
            yhat1      = DD1 + ZZ1*alphahat1 ;
            yhat2      = DD2 + ZZ2*alphahat1 ;
            yhat3      = DD3 + ZZ3*alphahat1 ;
            
            nu1        = YY(t,:) - yhat1';
            nu2        = YY(t,:) - yhat2';
            nu3        = YY(t,:) - yhat3';
            
            Ft1        = ZZ1*Phat1*ZZ1' + RR ;
            Ft2        = ZZ2*Phat1*ZZ2' + RR ;
            Ft3        = ZZ3*Phat1*ZZ3' + RR ;
           
            Ft1        = 0.5*(Ft1 + Ft1');
            Ft2        = 0.5*(Ft2 + Ft2');
            Ft3        = 0.5*(Ft3 + Ft3');
         
            % Compute likelihood
      
            loglht(i,1)      = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft1))) - 0.5*nu1/(Ft1)*nu1' );
            loglht(nx+i,1)   = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft2))) - 0.5*nu2/(Ft2)*nu2' );
            loglht(2*nx+i,1) = real( -0.5*size(YY,2)*log(2*pi)- 0.5*real(log(det(Ft3))) - 0.5*nu3/(Ft3)*nu3' );

            % Updating Step
     
            At(i,:)      = ( alphahat1 + Phat1*ZZ1'/(Ft1)*nu1' )';
            At(nx+i,:)   = ( alphahat1 + Phat1*ZZ2'/(Ft2)*nu2' )';
            At(2*nx+i,:) = ( alphahat1 + Phat1*ZZ3'/(Ft3)*nu3' )';
                  
            Pt(i,:)      = vec(Phat1 - Phat1*ZZ1'/(Ft1)*ZZ1*Phat1' )';
            Pt(nx+i,:)   = vec(Phat1 - Phat1*ZZ2'/(Ft2)*ZZ2*Phat1' )';
            Pt(2*nx+i,:) = vec(Phat1 - Phat1*ZZ3'/(Ft3)*ZZ3*Phat1' )';
                          
        end   % end for i
          
          %Ft            = Ft-yhat*yhat';

      % loglh computation
      loglhtmax = max(loglht);
      loglh     = loglh + loglhtmax + log(probx2'*exp(loglht-loglhtmax));
      
      % Updating mixture probabilities and the state probabilities

      probx2 = (probx2.*exp(loglht-loglhtmax)) / ( probx2' * exp(loglht-loglhtmax) );
      probs  = [sum(probx2(1:nx)) ; sum(probx2(nx+1:2*nx)) ;sum(probx2(2*nx+1:3*nx)) ];
      
              
      % Collapse the state vector      
      % Check for zero probability states
      for j=1:nx
      i2 = (j-1)*8+1;
      probx(j) =   probx2(i2) +probx2(i2+1) +probx2(i2+2);
            
      At1(j,:) =    (probx2(i2)/probx(j)) * At(i2,:)  + (probx2(i2+1)/probx(j)) * At(i2+1,:)...
                 +(probx2(i2+2)/probx(j)) * At(i2+2,:);
             
      Pt1(j,:) = vec( (probx2(i2)/probx(j)) * ( (reshape(Pt(i2,:),nstate,nstate)  +At(i2,:)'*At(i2,:)))...
                  + (probx2(i2+1)/probx(j)) * ( (reshape(Pt(i2+1,:),nstate,nstate)+At(i2+1,:)'*At(i2+1,:)))...
	  	          + (probx2(i2+2)/probx(j)) * ( (reshape(Pt(i2+2,:),nstate,nstate)+At(i2+2,:)'*At(i2+2,:)))...
                   - At1(j,:)'*At1(j,:))';
      end
      
      At_mat(t,:) = At1(1,:) ;
      probx = probx/sum(probx);
      probsT(t,:)=probs';
end

bubble = probsT(:,2).*(1+gx(end,:)*At_mat')' ;

At_mat2prob = [] ; 
    
%  % Backward simulation 
%    t = nobs; 
%          
%         probs1     = probsT(t,1);
%         probs2     = probsT(t,2);
%         
%         aux = rand(1,1);
%         
%         draws1     = aux <= probs1;
%         draws2     = aux <= probs1+probs2;
%         
%         if draws1+draws2==2;
%             ssT_draw(t,:)=[1 0 0];
%         elseif draws1+draws2==1;
%             ssT_draw(t,:)=[0 1 0];
%         elseif draws1+draws2==0;
%             ssT_draw(t,:)=[0 0 1];
%         end
%      
%      
%    loglh = log( ssT_draw(t,:)*probsT(t,:)');
% 
%    while t ~= 1;
%   
%       t = t-1;
%       probs1 = ssT_draw(t+1,:)*strans(:,1)*probsT(t,1)/( ssT_draw(t+1,:)*strans*probsT(t,:)' );
%       probs2 = ssT_draw(t+1,:)*strans(:,2)*probsT(t,2)/( ssT_draw(t+1,:)*strans*probsT(t,:)' );
%       
%         aux = rand(1,1);
%         
%         draws1     = aux <= probs1;
%         draws2     = aux <= probs1+probs2;
%            
%         if draws1+draws2==2;
%             ssT_draw(t,:)=[1 0 0];
%         elseif draws1+draws2==1;
%             ssT_draw(t,:)=[0 1 0];
%         elseif draws1+draws2==0;
%             ssT_draw(t,:)=[0 0 1];
%         end
%         
%         loglh = loglh + log( ssT_draw(t,:)*[probs1;probs2;(1-probs1-probs2)] );
%    
%    end   


