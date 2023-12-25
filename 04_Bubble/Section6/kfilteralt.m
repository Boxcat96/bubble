function [shatnew,signew,loglh, ypred]=kfilteralt(y,H,F,c,C,shat,sig,R,Q,n)

% shat = s{t|t-1}
% sig  = sigma{t|t-1}
% Hamilton page 381 and 385

res     = y - c - H*shat ;

temp    = H*sig*H' + R ;
% sigtemp = pinv(temp) ;
sigtemp = temp\eye(size(temp)) ;
K       = F*sig*H'*sigtemp ;
shatnew = C + F*shat + K*res ;
signew  = ((F - K*H)*sig*(F' - H'*K') + K*R*K') + Q ;

ypred   = c + H*(shat) ;

loglh=-.5*log((det(temp))) - .5*res'*sigtemp*res - .5*n*log(2*pi);
