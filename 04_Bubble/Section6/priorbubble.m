function logprior = priorgrowth(xx)

% Inputs:   AI: distribution description from pracel
%           MI: mean and variance from pracel
%           xx: vector of estimated coefficients
% Ouput:    log prior distribution
% Pablo A. Guerron, NCSU, 08 - Dec
% This uses Del Negro - Schorfheide Inverse Gamma definition
% get distributions , means, and variances

% load prsoe
 
limites ;
 
setup_bubble_prior ;

logprior = 0 ;
jj = 1;
for i = 1:size(II,1) ;
    if II(i,1) == 1 
        if AI(i,1) == 'I';
            % Next lines inverse gamma taken from Marco-Frank code
            a = MI(i,1); % enter first parameter in priors as it were mean
            b = MI(i,2); % second parameter as std
            logprior = logprior - gammaln(a) + log(exp(-b/xx(jj))*b^a*xx(jj)^(-1-a)) ; % Wolfram
            % logprior = logprior +log(2)-gammaln(b/2)+(b/2)*log(b*a^2/2)-((b+1)/2)*log(xx(jj)^2)-b*a^2/(2*xx(jj)^2);
            % logprior = logprior - log(xx(jj)) ; % Jeffrey's uninformative prior
            % Alternative inverse gamma
            % [m,n]   =   meaninvgammamod(MI(i,1),MI(i,2)) ;
            % logprior =  logprior - gammaln(m) + log(n^m*(1/xx(jj))^(m+1)*exp(-n/xx(jj))) ;         
            jj = jj + 1 ; 
        elseif AI(i,1) == 'N'
            % normal
            logprior = log(pdf('Normal',xx(jj),MI(i,1),MI(i,2))) + logprior ;
            jj = jj + 1 ;
        elseif AI(i,1) == 'B'
            % beta
            [p,q] = meanbeta(MI(i,1),MI(i,2)) ;
            logprior = log(pdf('Beta',xx(jj),p,q)) + logprior ;
            jj = jj + 1 ;
        elseif AI(i,1) == 'G'
            % gamma
            [a,b] = meangamma(MI(i,1),MI(i,2)) ;
            logprior = log(pdf('Gamma',xx(jj),a,b)) + logprior ;
            jj = jj + 1 ;
        elseif AI(i,1) == 'R'
            % gamma
            [a,b] = meangamma(MI(i,1),MI(i,2)) ;
            logprior = log(pdf('Gamma',xx(jj)-1,a,b)) + logprior ;
            jj = jj + 1 ;            
        elseif AI(i,1) == 'P'
            % sigw, sig, and sigm
            logprior = logprior + log(pdf('Normal',1/(xx(jj) - 1),MI(i,1),MI(i,2))) ; 
            jj = jj + 1 ;
        elseif AI(i,1) == 'U'
            logprior = logprior + log(pdf('Uniform',xx(jj),MI(i,1),MI(i,2))) ; 
            jj = jj + 1 ;            
        end
    end
end

function [p,q] = meanbeta(meanbet,varbet) 

p           =   (1 - meanbet)*meanbet^2/varbet - meanbet ;
q           =   p/meanbet - p  ;

function [a,b]  =   meangamma(meangam,vargam)

b           =   vargam/meangam ;
a           =   meangam/b ;

function [m,n]  =   meaninvgamma(invmean,invar) 

n           =   2*invmean^2/invar + 4 ;
m           =   invmean*(n - 2) ;

function [m,n]  =   meaninvgammamod(invmean,invar) 
% implement as in wiki
m           =   invmean^2/invar + 2 ;
n           =   invmean*(m - 1 ) ;