% anal_deriv_print2f.m
%anal_deriv_print2f(filename,fx,fxp,fy,fyp,f,ETASHOCK, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx);
% Writes symbolic derivatives of the equilibrium conditions of a DSGE model  to an m-file for numeric evaluation.
% inputs: 
%         filename: ' ...' name of the file to append to '_num_eval'
%        The output of the program anal_deriv.m
%        ETASHOCK: symbolic matrix scaling std dev of shocks (the evolution of the state  is: x_t+1 = hx x_t + ETASHOCK epsilon_t+1, where epsilon_t+1 has an identity var/cov matrix)
% output: m-file "filename_num_eval.m" 

% by Andrea Pescatori Jun 18 2008, modified Aug 08
function anal_deriv_print2f(filename,fx,fxp,fy,fyp,f, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx);

if nargin>8
    approx=2;
else
    approx =1;
end


fid=fopen([filename,'_num_eval.m'],'w');
fprintf(fid,'%% File name: %s_num_eval.m \n',filename);
fprintf(fid,'%% File generated by anal_deriv_print2f.m Date: %s\n\n',date);

S = char(fx(:));
S=['nfx=', S(8:end-1),';\n'];  
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(fxp(:));
S=['nfxp=', S(8:end-1),';\n']; 
S=regexprep(S, ',', '\r');
fprintf(fid,S);

S = char(fy(:));
S=['nfy=', S(8:end-1),';\n'];
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(fyp(:));
S=['nfyp=', S(8:end-1),';\n']; 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

S = char(f(:));
S=['nf=', S(8:end-1),';\n']; 
S=regexprep(S, ',', '\r'); 
fprintf(fid,S);

% S = char(ETASHOCK(:));
% S=['nETASHOCK=', S(8:end-1),';\n'];  
% S=regexprep(S, ',', '\r');
% fprintf(fid,S);

S2F = ['nfx=reshape(nfx,[', num2str(size(fx)),']);\n'];  % (1)
fprintf(fid,S2F);

S2F = ['nfxp=reshape(nfxp,[', num2str(size(fxp)),']);\n'];  % (1)
fprintf(fid,S2F);

S2F = ['nfy=reshape(nfy,[', num2str(size(fy)),']);\n'];  % (1)
fprintf(fid,S2F);

S2F = ['nfyp=reshape(nfyp,[', num2str(size(fyp)),']);\n'];  % (1)
fprintf(fid,S2F);

S2F = ['nf=reshape(nf,[', num2str(size(f)),']);\n\n'];  % (1)
fprintf(fid,S2F);

% S2F = ['nETASHOCK=reshape(nETASHOCK,[', num2str(size(ETASHOCK)),']);\n'];  % (1)
% fprintf(fid,S2F);


if approx==1,
    fclose(fid);
     return
end

S=char(fypyp,1); % vectorizing
S=['nfypyp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfypyp=reshape(nfypyp,[', num2str(size(fypyp)),']);\n'];   
fprintf(fid,SF);

S=char(fypy,1);
S=['nfypy=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfypy=reshape(nfypy,[', num2str(size(fypy)),']);\n'];  
fprintf(fid,SF);

S=char(fypxp,1);
S=['nfypxp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfypxp=reshape(nfypxp,[', num2str(size(fypxp)),']);\n'];  
fprintf(fid,SF);

S=char(fypx,1);
S=['nfypx=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfypx=reshape(nfypx,[', num2str(size(fypx)),']);\n'];  
fprintf(fid,SF);

S=char(fyyp,1);
S=['nfyyp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfyyp=reshape(nfyyp,[', num2str(size(fyyp)),']);\n'];  
fprintf(fid,SF);

S=char(fyy,1);
S=['nfyy=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfyy=reshape(nfyy,[', num2str(size(fyy)),']);\n'];  
fprintf(fid,SF);

S=char(fyxp,1);
S=['nfyxp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfyxp=reshape(nfyxp,[', num2str(size(fyxp)),']);\n'];  
fprintf(fid,SF);

S=char(fyx,1);
S=['nfyx=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfyx=reshape(nfyx,[', num2str(size(fyx)),']);\n'];  
fprintf(fid,SF);

S=char(fxpyp,1);
S=['nfxpyp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxpyp=reshape(nfxpyp,[', num2str(size(fxpyp)),']);\n'];  
fprintf(fid,SF);

S=char(fxpy,1);
S=['nfxpy=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxpy=reshape(nfxpy,[', num2str(size(fxpy)),']);\n'];  
fprintf(fid,SF);

S=char(fxpxp,1);
S=['nfxpxp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxpxp=reshape(nfxpxp,[', num2str(size(fxpxp)),']);\n'];  
fprintf(fid,SF);

S=char(fxpx,1);
S=['nfxpx=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxpx=reshape(nfxpx,[', num2str(size(fxpx)),']);\n'];  
fprintf(fid,SF);

S=char(fxyp,1);
S=['nfxyp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxyp=reshape(nfxyp,[', num2str(size(fxyp)),']);\n'];  
fprintf(fid,SF);

S=char(fxy,1);
S=['nfxy=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxy=reshape(nfxy,[', num2str(size(fxy)),']);\n'];  
fprintf(fid,SF);

S=char(fxxp,1);
S=['nfxxp=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxxp=reshape(nfxxp,[', num2str(size(fxxp)),']);\n'];  
fprintf(fid,SF);

S=char(fxx,1);
S=['nfxx=', S(8:end-1),';\n']; S=regexprep(S, ',', '\r');
fprintf(fid,S);
SF = ['nfxx=reshape(nfxx,[', num2str(size(fxx)),']);\n'];  
fprintf(fid,SF);

fclose(fid);
