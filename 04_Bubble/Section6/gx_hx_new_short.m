function [gx,hx,exitflag]=gx_hx_new_short(fy,fx,fyp,fxp,stake);
% global gdp
%[gx,hx,exitflag]=gx_hx(fy,fx,fyp,fxp,stake);
%computes the matrices gx and hx that define the first-order approximation 
%of the DSGE model. That is, if 
%E_t[f(yp,y,xp,x)=0, then the solution is of the form
%xp = h(x,sigma) + sigma * eta * ep
%y = g(x,sigma).
%The first-order approximations to the functions g and h around the point (x,sigma)=(xbar,0), where xbar=h(xbar,0), are:
%h(x,sigma) = xbar + hx (x-xbar) 
%and
%g(x,sigma) = ybar + gx * (x-xbar),
%where ybar=g(xbar,0). 
%%The exit flag takes the values 0 (no solution), 1 (unique solution), 2 (indeterminacy), or 3 (z11 is not invertible).
%Inputs: fy fyp fx fxp
%Outputs: gx hx exitflag
%Calls subfunction solab_flag.m 
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001, May 11 2006
A = [-fxp -fyp];
B = [fx fy];

if nargin<5
    stake=1;
end
% nancheck=sum(isnan([A(:) B(:)]));
% if nancheck==0
    [gx,hx,exitflag]=solab_flag(A,B,size(fx,2),stake);
% else
%     exitflag=4
%     gx=[];hx=[];
% end


%SOLAB_FLAG.M
% Function: solab_flag
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number NK
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% NK is the number of state variables.
%
% Outputs: the decision rule f, the law of motion p, and the exit flag exitflag.  If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
%The exit flag takes the values 0 (no solution), 1 (unique solution), 2 (indeterminacy), or 3 (z11 is not invertible).
% Calls: qzdiv.m
%
%This program is a slightly modified version of the program solab.m by Paul Klein. It allows for a more convenient handling of cases of no local equilibrium and multyple equilibria.
%
%(c) Last modified May 11, 2006

function [f,p,exitflag] = solab_flag(a,b,NK,stake);
global para1  PHI
if nargin<4
    stake = 1;
end


exitflag = 1;


[s,t,q,z] = qz(a,b);   % upper triangular factorization of the matrix pencil b-za

[s,t,q,z] = ordqz(s,t,q,z,'udo');
% 
% [s,t,q,z] = qzdiv(stake,s,t,q,z);   % reordering of generalized eigenvalues in ascending order

nk=sum(abs(diag(t))<stake*abs(diag(s))); %stake*

if nk>NK
%warning('The Equilibrium is Locally Indeterminate')

exitflag=2;

elseif nk<NK

%warning('No Local Equilibrium Exists')
exitflag = 0;
end

z21 = z(nk+1:end,1:nk);
z11 = z(1:nk,1:nk);

if rank(z11)<nk;
%	warning('Invertibility condition violated')
   exitflag = 3;
end

if exitflag==1; %if the solution is uniquely determinated, adn the invertibility conditions are satisfied
    
s11 = s(1:nk,1:nk);
t11 = t(1:nk,1:nk);

dyn = s11\t11;

z11i = z11\eye(nk);

f = real(z21*z11i);  % The real function takes away very small imaginary parts of the solution

p = real(z11*dyn*z11i);

else % if not, don't solve the model
    f=[];
    p=[];
end
    


