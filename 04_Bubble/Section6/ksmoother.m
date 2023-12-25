function [shatnew,signew,omega, ypred]=ksmoother(y,H,F,c,C,shat,sig,R,Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t)=c+H*s(t)+e(t)
% s(t)=C+F*s(t-1)+v(t)
% V(e(t))=R
% V(v(t))=Q
% Page 26 in state-space models by kim and nelson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(y);
omega=F*sig*F'+Q;
sigma=H*omega*H'+R;
isigm = sigma\eye(size(sigma)) ;
% k=omega*H'/sigma;
k = omega*H'*isigm;

res=y-c-H*(F*shat+C);
shatnew=C+F*shat+k*res;
signew=omega-k*H*omega;
% ypred = c + H*(F*shat + C) ;
ypred = c + H*(shatnew + C) ;

% loglh=-.5*log(det(sigma)) - .5*(res'/sigma)*res - .5*n*log(2*pi);
% loglh=-.5*log(det(sigma)) - .5*res'*isigm*res - .5*n*log(2*pi);