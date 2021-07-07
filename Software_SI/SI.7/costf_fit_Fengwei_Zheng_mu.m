%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Evaluation of the relationship mass-cell specific growth rate
% 
% we use the findings of Zheng etal., General quantitative relations linking 
%   cell growth and the cell cycle in Escherichia coli, Nature Microbiol.,
%   2020. DOI: ﻿10.1038/s41564-020-0717-x
% and
%  ﻿Fangwei Si etal. ,﻿Invariance of Initiation Mass and Predictability 
%   of Cell Size in Escherichia coli, ﻿Current Biology 27, 1278–1287, 2017.
%
% and data from Bremer and Dennis, Modulation of Chemical Composition and 
%      Other Parameters of the Cell by Growth Rate, 1996.
%
% We consider the relationships
%
%   1/(C+D) = mu/(alpha+beta*mu) from Zheng's work and
%   S = S_0 2^{(C+D)/t_d}, with t_d = log(2)/mu
%   m = density*S
% where:
% C= Time required to duplicate the chromosome
% D= Time ﻿time interval between termination of replication and completion
%     of cell division
% mu = specific growth rate
% t_d = duplication time
% S_0 = cell initiation volume 
% m_0 = 112e-15 (g/cell), initiation cell dry weight. Estimated from 
%       S_0=0.28e-18 m^3/cell (Fengwei, 2017), 
%       E coli cell density 1.1g/ml (Bionumbers)
%       E coli water content 70% (Bionumbers)
% S = cell volume
%
% Thus, we get the model: mc = m0*exp(alpha+beta*mu)
% that can be expreseed as:
%    log(mc) = log(m0) + alpha + beta*mu 
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
% This is the main MATLAB script. The cost function is implemented imn the
% script costf_fit_Fengwei_Zheng_mu.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This is the cost function MATLAB script. 
%  The main MATLAB script is implemented in Main_fit_Fengwei_Zheng_mu.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function J = costf_fit_Fengwei_Zheng_mu(x) 

% Global variables:
global cellMass;
global mu;
global m_0;

parameters = x; 
Prediction=zeros(size(mu))'; %Allocate space
e=zeros(size(mu))';  %Allocate space        
   
for k=1:length(mu)
  Prediction(k) =parameters(1)*exp(parameters(2)*mu(k)^parameters(3));
  %Prediction(k) = parameters(1)*exp(parameters(2)*mu(k));  
  e(k)=cellMass(k)-Prediction(k); 
end
   
J=0.5*e'*e; % Mean square prediction error
  





   
     
   
 
 
 