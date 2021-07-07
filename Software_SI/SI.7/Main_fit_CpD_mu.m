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
% This is the main MATLAB script. The cost function is implemented imn the
% script costf_fit_CpD_mu.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global variables:
global Measured;
global mu;

td=[24, 30, 40, 60 100]/60; %(hours)
mu =log(2)./td;

% Experimental data (hours):
C_exp = [42
    43
    45
    50
    67]/60;
    
D_exp = [23
    24
    25
    27
    30]/60;

Measured = 1./(C_exp+D_exp);

%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='costf_fit_CpD_mu';    %script cost function to optimize

problem.x_L=[0,0]; % Expected minimum values alpha, beta
problem.x_U=[2,2]; % Expected maximum values alpha, beta

opts.maxeval=10000; 
opts.maxtime=1e7;  
opts.local.solver='fminsearch';
opts.local.iterprint=1;

%========================= END OF PROBLEM SPECIFICATIONS =====================
% rand('state',0)
Results=MEIGO(problem,opts,'ESS');

parameters =Results.xbest; % Vector with best estimated values

%%%%%%% PLOT RESULTS  %%%%%%%%%%%

figure(1)
box_i = [0.3:0.1:1.8];
y=[];
for k=1:length(box_i)
  y(k,1) = box_i(k)/(parameters(1)+parameters(2)*box_i(k));  
end

   
plot(box_i,y,'b-'),ylabel('\mu/(\alpha+\beta*\mu)'), xlabel('Growth rate')
hold on
plot(mu,Measured,'rs','MarkerSize',10)
grid on
legend('Estimated \mu/(\alpha+\beta*\mu)', '1/(C+D)')
hold off




