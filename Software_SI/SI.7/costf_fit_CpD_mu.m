%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Evaluation of the relationship mass-cell specific growth rate
% 
% We use the findings of Zheng etal., General quantitative relations linking 
%   cell growth and the cell cycle in Escherichia coli, Nature Microbiol.,
%   2020. DOI: ﻿10.1038/s41564-020-0717-x
%
% and data from Bremer and Dennis, Modulation of Chemical Composition and 
%      Other Parameters of the Cell by Growth Rate, 1996.
%
% We consider the relationship
%
% 1/(C+D) = mu/(alpha+beta*mu)
% 
% where:
% C= Time required to duplicate the chromosome
% D= Time ﻿time interval between termination of replication and completion
%     of cell division
% mu = specific growth rate
% 
% Time units are in minutes
% 
% Parameters to estimate:
%					alpha
% 					beta
%
% To estimate the parameters we use the algorithm MEIGO: Egea JA, Henriques D, 
% Cokelaer T, Villaverde AF, MacNamara A, Danciu DP, Banga JR and Saez-Rodriguez J. 
% (2014) MEIGO: an open-source software suite based on metaheuristics for
% global optimization in systems biology and bioinformatics. 
% BMC Bioinformatics 15:136. 
% https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-136
%
%  This is the cost function MATLAB script. The main MATLAB script is implemented 
%  in the script main_fit_CpD_mu.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function J = costf_fit_CpD_mu(x) 

% Global variables:
global Measured;
global mu;

parameters = x; 
Prediction=zeros(size(mu))'; %Allocate space
e=zeros(size(mu))';  %Allocate space        
   
for k=1:length(mu)
  Prediction(k) = mu(k)/(parameters(1)+parameters(2)*mu(k));  
  e(k)=Measured(k)-Prediction(k); 
end
   
J=0.5*e'*e; % Mean square prediction error
  





   
     
   
 
 
 