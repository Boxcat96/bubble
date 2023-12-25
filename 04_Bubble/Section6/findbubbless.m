function ff = findbubbless(xsol,param) 

ifrac   = param.ifrac ;
alph    = param.alph ;
dEls    = param.dEls ;
phi     = param.phi ;
sigmab  = param.sigmab ;
sigmaf  = param.sigmaf ;
deprts  = param.deprts ;
labs    = param.labs ;
bet     = param.bet ;
rho     = param.rho ;
rs      = param.rs ;
utils   = param.utils ;

yft      = xsol(1) ;
ybt      = xsol(2) ;
uft      = xsol(3) ;
ubt      = xsol(4) ;
qft      = xsol(5) ;
qbt      = xsol(6) ;
cift     = xsol(7) ;
cibt     = xsol(8) ;
csft     = xsol(9) ;
csbt     = xsol(10) ;
labft    = xsol(11) ;
labbt    = xsol(12) ;
eta      = xsol(13) ;
delta    = xsol(14) ;
d1       = xsol(15) ;
mbt      = xsol(16) ;

lambdaft = ((qft) - 1)/(1 - phi*(qft));
lambdabt = ((qbt) - 1)/(1 - phi*(qbt));
rft      = alph*(yft)/(uft);
rbt      = alph*(ybt)/(ubt);
wft      = (1-alph)*(yft)/((1-ifrac)*(labft));
wbt      = (1-alph)*(ybt)/((1-ifrac)*(labbt));
cft      = ifrac*(cift) + (1-ifrac)*(csft);
cbt      = ifrac*(cibt) + (1-ifrac)*(csbt);
ivstft   = ifrac*((uft)*(rft) + phi*(qft)*(1 - (delta + d1*( (uft)-1 + dEls/2*( (uft)-1 )^2 ))))/(1-phi*(qft));
ivstbt   = ifrac*((ubt)*(rbt) + phi*(qbt)*(1- (delta + d1*( (ubt)-1 + dEls/2*( (ubt)-1 )^2 ))) + (mbt))/(1-phi*(qbt));
gft      = 1 - (delta + d1*( (uft)-1 + dEls/2*( (uft)-1 )^2 )) + (ivstft);
gbt      = 1 - (delta + d1*( (ubt)-1 + dEls/2*( (ubt)-1 )^2 )) + (ivstbt);

ff = zeros(16,1) ;

ff(1)  = -yft/ybt + (uft/ubt)^alph*(labft/labbt)^(1-alph) ;
ff(2)  = -(cift)^(-rho) + (csft)^(-rho)*((1-(labft))/(1-labs))^(eta*(1-rho)) ;
ff(3)  = -(cibt)^(-rho) + (csbt)^(-rho)*((1-(labbt))/(1-labs))^(eta*(1-rho)) ;
ff(4)  = -eta*csft/(1-labft) + wft ;
ff(5)  = -eta*csbt/(1-labbt) + wbt ;
ff(6)  = rft - d1*(1 + dEls*(uft - 1))*qft + ifrac*lambdaft*(rft - phi*qft*d1*(1 + dEls*(uft - 1))) ;
ff(7)  = rbt - d1*(1 + dEls*(ubt - 1))*qbt + ifrac*lambdabt*(rbt - phi*qbt*d1*(1 + dEls*(ubt - 1))) ;
ff(8)  = -qft + (1 - sigmaf)*bet*(1/gft)^rho*(uft*rft + (1 - (delta + d1*(uft - 1 + dEls/2*(uft - 1)^2)))*qft ...
                + ifrac*lambdaft*(uft*rft + phi*qft*(1 - (delta + d1*(uft - 1 + dEls/2*(uft - 1)^2))))) ...
                + sigmaf*bet*((cift/cibt)/gft)^rho*(ubt*rbt + (1 - (delta + d1*(ubt - 1 + dEls/2*(ubt - 1)^2)))*qbt ...
                + ifrac*lambdabt*(ubt*rbt + phi*qbt*(1 - (delta + d1*(ubt - 1 + dEls/2*(ubt - 1)^2))))) ;             
ff(9)  = -qbt + (1 - sigmab)*bet*(1/gbt)^rho ...
            *(ubt*rbt + (1 - (delta + d1*(ubt - 1 + dEls/2*(ubt - 1)^2)))*qbt ...
                + ifrac*lambdabt*(ubt*rbt + phi*qbt*(1 - (delta + d1*(ubt - 1 + dEls/2*(ubt - 1)^2))))) ...
                 + sigmab*bet*((cibt/cift)/gbt)^rho ...
            *(uft*rft + (1 - (delta + d1*(uft - 1 + dEls/2*(uft - 1)^2)))*qft ...
                + ifrac*lambdaft*(uft*rft + phi*qft*(1 - (delta + d1*(uft - 1 + dEls/2*(uft - 1)^2))))) ;
ff(10) = -(yft) + (cft) + (ivstft);
ff(11) = -(ybt) + (cbt) + (ivstbt);
ff(12) = sigmab/(sigmab+sigmaf)*labft + sigmaf/(sigmab+sigmaf)*labbt - labs ;
ff(13) = sigmab/(sigmab+sigmaf)*uft   + sigmaf/(sigmab+sigmaf)*ubt - utils ;
ff(14) = sigmab/(sigmab+sigmaf)*rft   + sigmaf/(sigmab+sigmaf)*rbt -rs ;
ff(15) = sigmab/(sigmab+sigmaf)*((delta + d1*(uft - 1 + dEls/2*(uft - 1)^2))) + sigmaf/(sigmab+sigmaf)*(delta + d1*(ubt - 1 + dEls/2*(ubt - 1)^2)) - deprts;
ff(16) = -1 + (1 - sigmab)*bet*(1/gbt)^rho*( 1 + ifrac*lambdabt)*gbt;