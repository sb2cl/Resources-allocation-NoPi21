%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               
% Evaluation of the relationship protein & cell dry weight-cell specific growth rate
% 
% we use the findings of Zheng etal., General quantitative relations linking 
%   cell growth and the cell cycle in Escherichia coli, Nature Microbiol.,
%   2020. DOI: ﻿10.1038/s41564-020-0717-x
%
% and data from Bremer and Dennis, Modulation of Chemical Composition and 
%      Other Parameters of the Cell by Growth Rate, 1996.
%
% We considered the relationships:
%  
%        dm_p/dmu = beta*m_p
% so that:
%        m_p = m_p0*exp(beta*mu)
%
% and 
%       dm_c/dmu = gamma*m_c
% so that:
%        m_c = m_c0*exp(gamma*mu)
% where:
%       mc = cell DW
%       mp = cell protein DW
%       mu = specific growth rate
%
%  Results obtained using Main_fit_FZ_mu_mp.m and Main_fit_Fengwei_Zheng_mu.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Experimental data (mins):
td_exp=[24, 30, 40, 60 100];
mu =log(2)./td_exp;

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

% Our best estimated parameters using MATLAB script main_fit_Fengwei_Zheng_mu.m
% for mc = m0*exp(gamma*mu):
     m0 = 123.022e-15; %(g)
     gamma = 68.5547; % (mins)   

% Our best estimated parameters using MATLAB script main_fit_Fengwei_Zheng_mu.m
% for mp = mp0*exp(beta*mu):
     mp0 = 77.3748e-15; % (g)
     beta = 61.7813; % (mins)   

% Notice mp0/m0 = 77.3748e-15/123.022e-15 = 63%    

box_mp = [0.5*min(proteinMass):5e-15:1.25*max(proteinMass)];
box_mc = [0.5*min(cellMass):5e-15:1.25*max(cellMass)];
box_mu = [0*min(mu):0.001:1.1*max(mu)];

f1=figure(1)
subplot(121)   
plot(box_mu,123.022e-15*exp(68.5547*box_mu),'b-','LineWidth',3)
xlabel('Growth rate (min^{-1})'), ylabel('Cell DW')
hold on
plot(mu,cellMass,'sk','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated m_c', 'Measured m_c','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off
subplot(122)   
plot(box_mu,77.3748e-15*exp(61.7813*box_mu),'b-','LineWidth',3),xlabel('Growth rate (min^{-1})'), ylabel('Protein DW')
hold on
plot(mu,proteinMass,'sk','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated m_p', 'Measured m_p','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off

% %f1.WindowState = 'maximized';
% %exportgraphics(f1,'../Texte/images/﻿f_mu_mc_mp.png','Resolution',300)

f2=figure(2) 
plot(mu,(proteinMass./cellMass),'rs','MarkerSize',10)
xlabel('Growth rate (min^{-1})'), ylabel('Protein/Cell DW')
grid on
hold on
plot(box_mu,(77.3748*exp(61.7813*box_mu))./(123.022*exp(68.5547*box_mu)),'b-','LineWidth',3)
ax = gca;
ax.FontSize = 14;
hold off

