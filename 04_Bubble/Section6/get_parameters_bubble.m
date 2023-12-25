% Parameters for Growth Model

% 1. Set parameter values

param.A     = 1 ;
param.eta   = 0.9 ;    
param.grows = 0.67 ; 
param.ifrac = 0.06 ;
param.alph  = 0.33 ;    
param.d1    = 0.03 ;   
param.dEls  = 0.33 ;
param.phi   = 0.19 ;
aax = 1/(1-1/147) ; % duration pre-Great Moderation; bubble
bbx = 1/(1-1/132) ;
p11 = (2-aax-(2-bbx)*(1-aax))/(1-(1-bbx)*(1-aax)) ;
p22 = 2 - bbx - (1-bbx)*p11 ;
% good results with 0.025
param.sigmab= 0.025 ; % 0.025 ; 0.05 ; 1 - p11 ; 0.01 ; 
param.sigmaf= 0.025 ; % 0.025 ; 0.05 ; 1 - p22 ; 0.01 ;
param.deprts= 0.025 ;
param.labs  = 0.25 ;
param.bet   = 0.99 ;
param.rho   = 1 ;
param.rs    = 0.05 ;

param.rhoa  = 0.95 ;
param.rhob  = 0.25 ;

param.siga = 0.005 ;
param.sigb = 0.005 ;

param.utils = 1.0 ;
param.utils = 1.0 ;

param.DD11 = 0.5 ;
param.DD12 = -0.05 ;
param.DD21 = 0.75 ;
param.DD22 = 0.05 ;
param.DD31 = 1.00 ;
param.DD32 = 0.025 ;

% compute some preliminaries                  
A       = param.A ;
eta     = param.eta ;
grows   = param.grows ;
ifrac   = param.ifrac ;
alph    = param.alph ;
d1      = param.d1 ;
dEls    = param.dEls ;
phi     = param.phi ;
sigmab  = param.sigmab ;
sigmaf  = param.sigmaf ;
deprts  = param.deprts ;
labs    = param.labs ;
bet     = param.bet ;
rho     = param.rho ;
yFs     = param.rs/alph ;
rhoa    = param.rhoa ;
rhob    = param.rhob ;
rs      = param.rs ;
utils   = param.utils ;
siga    = param.siga ;
sigb    = param.sigb ;
DD11    = param.DD11 ;
DD12    = param.DD12 ;
DD21    = param.DD21 ;
DD22    = param.DD22 ;
DD31    = param.DD31 ;
DD32    = param.DD32 ;