% function compute_steady_state_growth

% Compute steady state Growth model
% Pablo A. Guerron-Quintana
% Boston College
% Feb 08, 2017
% Revision Nov, 2018

% get_parameters_bubble ;


param1.A     = A ;
param1.eta   = eta ;    
param1.grows = grows ; 
param1.ifrac = ifrac ;
param1.alph  = alph ;    
param1.d1    = d1 ;   
param1.dEls  = dEls ;
param1.phi   = phi ;
param1.sigmab= sigmab ; 
param1.sigmaf= sigmaf ;
param1.deprts= deprts ;
param1.labs  = labs ;
param1.bet   = bet ;
param1.rho   = rho ;
param1.rs    = rs ;
param1.utils = utils ;

xsol = [0.11 0.13 0.86 1.13 1.63 1.30 0.08 0.09 0.08 0.09 0.24 0.25 1 0.025 1 1] ;

% benchmark
[xsol,fval,exitflag]    = fsolve(@(xx) findbubbless(xx,param1),xsol,optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000)) ;

if exitflag ~= 1 || xsol(16) < 0
    disp('No steady state')
    exitflags = 4 ;
    return
end

% findbubbless(xsol,param) ;

yft      = xsol(1) ;
yftp     = yft ;
ybt      = xsol(2) ;
ybtp     = ybt ;
uft      = xsol(3) ;
uftp     = uft ;
ubt      = xsol(4) ;
ubtp     = ubt ;
qft      = xsol(5) ;
qftp     = qft ;
qbt      = xsol(6) ;
qbtp     = qbt ;
cift     = xsol(7) ;
ciftp    = cift ;
cibt     = xsol(8) ;
cibtp    = cibt ;
csft     = xsol(9) ;
csftp    = csft ;
csbt     = xsol(10) ;
csbtp    = csbt ;
labft    = xsol(11) ;
labftp   = labft ;
labbt    = xsol(12) ;
labbtp   = labbt ;
eta      = xsol(13) ;
delta    = xsol(14) ;
d1       = xsol(15) ;
mbt      = xsol(16) ;
mbtp     = mbt ;

A        = (yft)/((uft^alph)*((1-ifrac)*labft)^(1-alph)) ;

lambdaft    = ((qft) - 1)/(1 - phi*(qft));
lambdaftp   = lambdaft ;
lambdabt    = ((qbt) - 1)/(1 - phi*(qbt));
lambdabtp   = lambdabt ;
rft         = alph*(yft)/(uft);
rftp        = rft ;
rbt         = alph*(ybt)/(ubt);
rbtp     = rbt ;
wft      = (1-alph)*(yft)/((1-ifrac)*(labft));
wftp     = wft ;
wbt      = (1-alph)*(ybt)/((1-ifrac)*(labbt));
wbtp     = wbt ;
cft      = ifrac*(cift) + (1-ifrac)*(csft);
cftp     = cft ;
cbt      = ifrac*(cibt) + (1-ifrac)*(csbt);
cbtp     = cbt ;
ivstft   = ifrac*((uft)*(rft) + phi*(qft)*(1 - (delta + d1*( (uft)-1 + dEls/2*( (uft)-1 )^2 ))))/(1-phi*(qft));
ivstftp  = ivstft ;
ivstbt   = ifrac*((ubt)*(rbt) + phi*(qbt)*(1- (delta + d1*( (ubt)-1 + dEls/2*( (ubt)-1 )^2 ))) + (mbt))/(1-phi*(qbt));
ivstbtp  = ivstbt ;
gft      = 1 - (delta + d1*( (uft)-1 + dEls/2*( (uft)-1 )^2 )) + (ivstft);
gftp     = gft ;
gbt      = 1 - (delta + d1*( (ubt)-1 + dEls/2*( (ubt)-1 )^2 )) + (ivstbt);
gbtp     = gbt ;

priskft  = (1 - sigmaf)*bet*(1/gft)^rho + sigmaf*bet*((cift/cibtp)/gft)^rho ;
priskftp = priskft ;
priskbt  = (1 - sigmab)*bet*(1/gbt)^rho + sigmab*bet*((cibt/ciftp)/gbt)^rho ;
priskbtp = priskbt ;

spreadft = 1/priskft - 1/qft ;
spreadftp= spreadft ;
spreadbt = 1/priskbt - 1/qbt ;
spreadbtp= spreadbt ;

tvft = ifrac*qft*phi*((1 - (delta + d1*( (uft)-1 + dEls/2*( (uft)-1 )^2 ))) + ivstft)/yft ;

tvbt = ifrac*(qbt*phi*((1 - (delta + d1*( (ubt)-1 + dEls/2*( (ubt)-1 )^2 ))) + ivstbt) + mbt)/ybt ;

turnft = phi*qft*gft/(yft) ;

turnbt = (phi*qbt*gbt + mbt)/(ybt) ;

ss_to_logs ;