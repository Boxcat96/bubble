function [var_y, var_x] = growth_var_decomp(xx,par,maxop)

if maxop == 1
    xx  =   map(xx) ;
end

if imag(xx(:)')*imag(xx(:)) ~= 0
    disp('error') ;
    return
end

bubbleundo ;

compute_steady_state_growth;

bubblegrow_num_eval ;

ne = 2 ; % four shocks
nstates = size(nfx,2) ;
npsvar  = size(nfy,2) ;

if abs(nf(:)'*nf(:))> 1e-11 || isinf(([nfx(:);nfxp(:);nfy(:);nfyp(:)]'* ...
        [nfx(:);nfxp(:);nfy(:);nfyp(:)])) == 1 ...
        || isnan(([nfx(:);nfxp(:);nfy(:);nfyp(:);]'* ...
        [nfx(:);nfxp(:);nfy(:);nfyp(:)])) == 1 ; 
    gx = zeros(npsvar,nstates)  ;
    hx = zeros(nstates,nstates) ;
    exitflags = 4 ; % error in numerical evaluation
    return
end

%First-order approximation
[gx,hx,~] = gx_hx_new_short(nfy,nfx,nfyp,nfxp) ;

% number of structural shocks
ETAMATRIX = zeros(size(hx,1)); 

for ixx = 1:(ne+1)
    
    ETAMATRIX(end-1,end-1)  = (ixx==1 || ixx == 3)*siga^2 ;    
    ETAMATRIX(end,end)      = (ixx==2 || ixx == 3)*sigb^2 ;    

    [var_y(:,:,ixx),var_x(:,:,ixx)] = mom(gx,hx,ETAMATRIX,0);

end