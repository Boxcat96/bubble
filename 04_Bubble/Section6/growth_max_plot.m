function lik = growth_max_plot(xx,par,zzd,~,maxop)
% % Solve model and compute likelihood
% % PaGq
% % Federal Reserve Bank of Phildelphia
% % Feb 1, 2013
% % Revision Dec 2014
% 

if maxop == 1
    xx  =   map(xx) ;
end
 
if imag(xx(:)')*imag(xx(:)) ~= 0
    lik = 10000000000 ;
    return
end
 
% Unwrap parameters
bubbleundo ;

% For counterfactuals uncomment next lines and save gx hx.
% sigmab = 0.00000001;
% sigmaf = 0.00000001 ;

% Compute steady state
exitflags = 0 ;
compute_steady_state_growth

if exitflags == 4
    lik = 100000000 ;
    return
end

% numerical evaluation of model
bubblegrow_num_eval ; 

ne = 2 ; % two shocks
nstates = size(nfx,2) ;
npsvar  = size(nfy,2) ;

if abs(nf(:)'*nf(:))> 1e-11 || isinf(([nfx(:);nfxp(:);nfy(:);nfyp(:)]'* ...
        [nfx(:);nfxp(:);nfy(:);nfyp(:)])) == 1 ...
        || isnan(([nfx(:);nfxp(:);nfy(:);nfyp(:);]'* ...
        [nfx(:);nfxp(:);nfy(:);nfyp(:)])) == 1 ; 
    gx = zeros(npsvar,nstates)  ;
    hx = zeros(nstates,nstates) ;
    exitflags = 4 ; % error in numerical evaluation
    lik = 100000000 ;
    % disp('error nf in ss')
    return
end

% number of structural shocks
ETAMATRIX = zeros(size(nfx,2),ne); 
ETAMATRIX(end-1,1)  = siga ;
ETAMATRIX(end,2)    = sigb ;

% solve gx hx
[gx,hx,exitflags] = gx_hx_new_short(nfy,nfx,nfyp,nfxp) ;

if exitflags ~= 1
    lik = 100000000 ;
    return
end

% plot IRFs
% plot_irfs_growth ;

% evaluate likelihood
% [lik,~]   = linear_growth(gx,hx,nstates,npsvar,ne,ETAMATRIX,t83,zzd) ;

phi = [1-sigmaf 1-sigmab] ;
dds = [DD11 DD12 DD21 DD22] ;

plotbubbleregime 