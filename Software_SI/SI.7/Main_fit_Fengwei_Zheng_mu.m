clear mex;
clear all;
close all;

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

% Global variables:
global cellMass;
global mu;
global m_0;

td=[24, 30, 40, 60 100]; %(hours)
mu =log(2)./td;

% Experimental data:
% Protein dry mass (g)
proteinMass = [450
    340
    234
    156
    100]*1e-15;

% cell DW (g)
cellMass = [
    865
     641
     433
     258
    148]*1e-15;

m_0 = 112e-15; %(g/cell), initiation cell dry weight.  

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='costf_fit_Fengwei_Zheng_mu';    %script cost function to optimize

 problem.x_L=[1e-16, 10,0.1]; % Expected minimum values alpha, beta
 problem.x_U=[2e-13, 50,1.0]; % Expected maximum values alpha, beta
% problem.x_L=[1e-14, 50]; % Expected minimum values alpha, beta
% problem.x_U=[2e-13, 75]; % Expected maximum values alpha, beta


opts.maxeval=10000; 
opts.maxtime=1e7;  
opts.local.solver='fminsearch';
opts.local.iterprint=1;

%========================= END OF PROBLEM SPECIFICATIONS =====================
% rand('state',0)
Results=MEIGO(problem,opts,'ESS');

parameters =Results.xbest; % Vector with best estimated values

%%%%%%% PLOT RESULTS  %%%%%%%%%%%

box_mu = [0*min(mu):0.0001:1.1*max(mu)];

figure(1)
estimated_mc=[];
for k=1:length(box_mu)
  estimated_mc(k,1) =parameters(1)*exp(parameters(2)*box_mu(k)^parameters(3));  
  %estimated_mc(k,1) = parameters(1)*exp(parameters(2)*box_mu(k));  

end   
plot(box_mu,estimated_mc,'b-','LineWidth',3),ylabel('m_c (DW)'), xlabel('Growth rate')
hold on
plot(mu,cellMass,'rs','MarkerSize',10)
grid on
legend('Estimated cell DW', 'Measured cell DW')
ax = gca;
ax.FontSize = 13;
hold off

%%% RESULTS OBTAINED MODEL TWO PARAMETERS:
% Best solution value		2.96147e-27
% Decision vector
% 	1.23022e-13
% 	1.14258*60=68.5547 (mins)

%%% RESULTS OBTAINED MODEL THREE PARAMETERS:
% Best solution value		2.50967e-29
% Decision vector
% 	1.9535e-16
% 	15.2264
% 	0.16788



