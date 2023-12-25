%this routine allows you to determine which parameters are to be estimated,
%and which are fixed. 
%II(ii,jj), ii = 1 means estimate this parameter, 0 means keep it fixed
%the value of jj matters if the associated parameter is to be estimated
%jj = 1, restrict to unit interval
%jj = 2, restrict to unit circle
%jj = 3, restrict to be positive
%jj = 4, restrict bigger than unity
%jj = 5, no restriction

%the code takes into account the setting for criterion in main.m (this
%determines whether you want all shocks in the criterion, or just a subset)

% check setup_bubble_prior if adding more estimated parameters

II(1,1) = 0;II(1,2)=3;  %A  
II(2,1) = 0;II(2,2)=1;  %eta  
II(3,1) = 0;II(3,2)=3;  %grows 
II(4,1) = 0;II(4,2)=1;  %ifrac
II(5,1) = 0;II(5,2)=1;  %alph
II(6,1) = 0;II(6,2)=3;  %d1
II(7,1) = 0;II(7,2)=3;  %dEls
II(8,1) = 0;II(8,2)=1;  %phi
II(9,1) = 0;II(9,2)=1;  %sigmab    
II(10,1) = 0;II(10,2)=1;%sigmaf
II(11,1) = 0;II(11,2)=1;%deprts
II(12,1) = 0;II(12,2)=3;%labs  
II(13,1) = 0;II(13,2)=1;%bet
II(14,1) = 0;II(14,2)=3;%rho
II(15,1) = 0;II(15,2)=3;%rs
II(16,1) = 0;II(16,2)=1;%utils
II(17,1) = 1;II(17,2)=1;%rhoa
II(18,1) = 1;II(18,2)=1;%rhob
II(19,1) = 1;II(19,2)=3;%siga
II(20,1) = 1;II(20,2)=3;%sigB
II(21,1) = 1;II(21,2)=5;%DD11 
II(22,1) = 1;II(22,2)=5;%DD12
II(23,1) = 1;II(23,2)=5;%DD21 --
II(24,1) = 1;II(24,2)=5;%DD22 --
II(25,1) = 0;II(25,2)=5;%DD31 --
II(26,1) = 0;II(26,2)=5;%DD32  --

name(1,:)   =   'A        ';
name(2,:)   =   'eta      ';
name(3,:)   =   'grows    ';
name(4,:)   =   'ifrac    ';
name(5,:)   =   'alph     ';
name(6,:)   =   'd1       ';
name(7,:)   =   'dEls     ';
name(8,:)   =   'phi      ';
name(9,:)   =   'sigmab   ';
name(10,:)  =   'sigmaf   ';
name(11,:)  =   'deprts   ';
name(12,:)  =   'labs     ';
name(13,:)  =   'bet      ';
name(14,:)  =   'rho      ';
name(15,:)  =   'rs       ';
name(16,:)  =   'utils    ';
name(17,:)  =   'rhoa     ';
name(18,:)  =   'rhob     ';
name(19,:)  =   'siga     ';
name(20,:)  =   'sigb     ';
name(21,:)  =   'DD11     ';
name(22,:)  =   'DD12     ';
name(23,:)  =   'DD21     ';
name(24,:)  =   'DD22     ';
name(25,:)  =   'DD31     ';
name(26,:)  =   'DD32     ';

%jj = 1, restrict to unit interval
%jj = 2, restrict to unit circle
%jj = 3, restrict to be positive
%jj = 4, restrict to be larger than one
%jj = 5, no restriction

inunitcircle=[];
nonneg=[];
unrestr=[];
inzeroone=[];
morethanoneohone=[];

ix=0;
[n,m]=size(II);
for ii = 1:n
    if II(ii,1) == 1
        ix=ix+1;
        if II(ii,2)== 1
            inzeroone=[inzeroone,ix];
        elseif II(ii,2) == 2
            inunitcircle=[inunitcircle, ix];            
        elseif II(ii,2) == 3
            nonneg=[nonneg, ix];
        elseif II(ii,2) == 4
            morethanoneohone=[morethanoneohone, ix];   
        elseif II(ii,2) == 5
            unrestr=[unrestr, ix];
        end
    end
end
        
save restr inzeroone inunitcircle nonneg morethanoneohone unrestr 

fid=fopen('bubbledo.m','w');
ix=0;
for ii = 1:n
    if II(ii,1) == 1
        ix=ix+1;
        dd=['xx(%d)  =  ',name(ii,:),';\n'];
        fprintf(fid,dd,ix);
    end
end
ix=0;
for ii = 1:n
    if II(ii,1) == 0
        ix=ix+1;
        dd=['par(%d)  =  ',name(ii,:),';\n'];
        fprintf(fid,dd,ix);
    end
end
fclose(fid);
fid=fopen('bubbleundo.m','w');
ix=0;
for ii = 1:n
    if II(ii,1) == 1
        ix=ix+1;
        dd=[name(ii,:),'  =  xx(%d) ;\n'];
        fprintf(fid,dd,ix);
    end
end
ix=0;
for ii = 1:n
    if II(ii,1) == 0
        ix=ix+1;
        dd=[name(ii,:),'  =  par(%d) ;\n'];
        fprintf(fid,dd,ix);
    end
end
fclose(fid);

fid = fopen('limites.m','w') ;
for ii = 1:n
   dd = ['II(%d,1) = ', int2str(II(ii,1)),'; '] ;
   fprintf(fid,dd,ii) ;
   dd = ['II(%d,2) = ', int2str(II(ii,2)),';\n'] ;
   fprintf(fid,dd,ii) ;   
end
fclose(fid) ;

setup_bubble_prior ;