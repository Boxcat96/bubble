% prior distributions for SOE model
% AI: distribution type
% MI(:,1) : mean 
% MI(:,2) : variance
% P special normal distribution for sig, sigm, sigw
% V no distribution
% I inverse gamma
% G gamma
% R gamma for growth rates 
% B beta
% U uniform between MI(1,1) and MI(1,2)
% Pablo A. Guerron
% Philadelphi Fed, August 2013 
% Modified August 2018

% distribution      prior mean and variance except for normals
% Warning: Only priors 17 - 26 are valid here! 

AI(1,1) = 'B';   MI(1,1) = 0.3;     MI(1,2) = 0.05^2 ;  %alph   
AI(2,1) = 'B';   MI(2,1) = 0.98;    MI(2,2) = 0.005^2 ; %bet  
AI(3,1) = 'B';   MI(3,1) = 0.5;     MI(3,2) = 0.2^2 ;   %capA     
AI(4,1) = 'B';   MI(4,1) = 0.2;     MI(4,2) = 0.2^2 ;   %sgs  
AI(5,1) = 'B';   MI(5,1) = 0.2;     MI(5,2) = 0.2^2  ;  %d0 
AI(6,1) = 'B';   MI(6,1) = 0.2;     MI(6,2) = 0.2^2 ;   %d1
AI(7,1) = 'G';   MI(7,1) = 0.1  ;   MI(7,2) = 0.025^2 ;   %d2
AI(8,1) = 'B';   MI(8,1) = 0.08;    MI(8,2) = 0.2^2 ;   %sy
AI(9,1) = 'B';   MI(9,1) = 0.8;     MI(9,2) = 0.2^2 ;   %utils    
AI(10,1) = 'G';  MI(10,1) = 1.6;    MI(10,2) = 0.2^2 ;  %omega --
AI(11,1) = 'B';  MI(11,1) = 0.33;   MI(11,2) = 0.20^2 ; %psia
AI(12,1) = 'G';  MI(12,1) = 1.1;    MI(12,2) = 0.2^2 ;  %sige
AI(13,1) = 'G';  MI(13,1) = 1.2;    MI(13,2) = 0.2^2 ;  %nu  
AI(14,1) = 'G';  MI(14,1) = 0.15;   MI(14,2) = 0.025^2 ; %theta 
AI(15,1) = 'G';  MI(15,1) = 0.15;   MI(15,2) = 0.025^2 ;  %phiks 
AI(16,1) = 'G';  MI(16,1) = 1.6;    MI(16,2) = 0.2^2 ;  %curlyeps
AI(17,1) = 'B';  MI(17,1) = 0.75;    MI(17,2) = 0.05^2 ;    %kappa     
AI(18,1) = 'B';  MI(18,1) = 0.15;    MI(18,2) = 0.05^2 ;  %grows   
AI(19,1) = 'I';  MI(19,1) = 6;      MI(19,2) = 1 ;      %deltan
AI(20,1) = 'I';  MI(20,1) = 6;      MI(20,2) = 1 ;      %chi        
AI(21,1) = 'N';  MI(21,1) = 0.5;    MI(21,2) = 0.1 ;  %eta  
AI(22,1) = 'N';  MI(22,1) = 0;      MI(22,2) = 1 ;  %ppsi
AI(23,1) = 'N';  MI(23,1) = 0.75;   MI(23,2) = 0.1 ;  %invsh
AI(24,1) = 'N';  MI(24,1) = 10 ;    MI(24,2) = 1 ;  %rdsh
AI(25,1) = 'N';  MI(25,1) = 1.00;   MI(25,2) = 0.01 ;  %lsh
AI(26,1) = 'N';  MI(26,1) = 0.025;   MI(26,2) = 0.01 ;  %taup
% AI(27,1) = 'B';  MI(27,1) = 0.40;   MI(27,2) = 0.2^2 ;  %taul
% AI(28,1) = 'B';  MI(28,1) = 0.5;    MI(28,2) = 0.2^2 ;  %rhosw 
% AI(29,1) = 'B';  MI(29,1) = 0.5;    MI(29,2) = 0.2^2 ;  %rhou  
% AI(30,1) = 'B';  MI(30,1) = 0.5;    MI(30,2) = 0.2^2 ;  %rhomei  
% AI(31,1) = 'B';  MI(31,1) = 0.5;    MI(31,2) = 0.2^2 ;  %rhoa  
% AI(32,1) = 'B';  MI(32,1) = 0.5;    MI(32,2) = 0.2^2 ;  %rhophik  
% AI(33,1) = 'B';  MI(33,1) = 0.5;    MI(33,2) = 0.2^2 ;  %rhog  
% AI(34,1) = 'B';  MI(34,1) = 0.5;    MI(34,2) = 0.2^2 ;  %rhopsia  
% AI(35,1) = 'B';  MI(35,1) = 0.5;    MI(35,2) = 0.2^2 ;  %rhochi  
% AI(36,1) = 'B';  MI(36,1) = 0.5;    MI(36,2) = 0.2^2 ;  %rhobet  
% AI(37,1) = 'I';  MI(37,1) = 6;      MI(37,2) = 1 ;      %sigmasw
% AI(38,1) = 'U';  MI(38,1) = 0;      MI(38,2) = 0.01 ;   %sigmau
% AI(39,1) = 'U';  MI(39,1) = 0;      MI(39,2) = 0.03 ;   %sigmam
% AI(40,1) = 'I';  MI(40,1) = 6;      MI(40,2) = 1 ;      %sigmaa
% AI(41,1) = 'I';  MI(41,1) = 6;      MI(41,2) = 1 ;      %sigmapk      
% AI(42,1) = 'I';  MI(42,1) = 6;      MI(42,2) = 1 ;      %sigmag
% AI(43,1) = 'I';  MI(43,1) = 6;      MI(43,2) = 1 ;      %sigmaps
% AI(44,1) = 'I';  MI(44,1) = 6;      MI(44,2) = 1 ;      %sigmac
% AI(45,1) = 'I';  MI(45,1) = 6;      MI(45,2) = 1 ;      %sigmab