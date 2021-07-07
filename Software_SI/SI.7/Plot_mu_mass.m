%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% We consider Zheng's etal. relationship:
%     mu/(alpha+beta*mu)=1/(C+D)
%
%     mc/(c0+c1*mc) = mu
%     mp/(d0+d1*mp) = mu
% where:
%       mc = cell DW
%       mp = cell protein DW
%       mu = specific growth rate
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% Experimental data (mins):
td_exp=[24, 30, 40, 60 100];
mu =log(2)./td_exp;

% Time required to duplicate the chromosome, C (mins):
C_exp = [42
    43
    45
    50
    67];

% Time ﻿time interval between termination of replication and completion
%  of cell division, D (mins):
D_exp = [23
    24
    25
    27
    30];

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

% Fit parameters from Zheng etal for mu/(alpha+beta*mu)=1/(C+D):
alpha_Z=0.28;
beta_Z_h= 0.99; %(hours)
beta_Z_m= beta_Z_h*60; %(mins)

% Our best estimated parameters using MATLAB script main_fit_CpD_mu.m
% for mu/(alpha+beta*mu)=1/(C+D)
	alpha= 0.277079;
	beta_h= 0.912718; %(hours)
    beta_m= beta_h*60; %(mins)
 
% Our best estimated parameters using MATLAB script main_fit_mc_mu.m
% for mp/(c0+c1*mp)=mu:
    c0 = 1.2364e-11; %(g*min)
    c1 = 6.959895861;    % (min)
       
   % Our best estimated parameters using MATLAB script main_fit_mc_mu.m
% for mc/(d0+d1*mc) = mu:
    d0 = 1.9641e-11; %(g*min)
    d1 = 12.1265289;    % (min)
    
    % Our best estimated parameters using MATLAB script main_fit_Fengwei_Zheng_mu.m
% for mc = m0*exp(gamma*mu):
    m0 = 123.022e-15; %(g)
    gamma = 1.14258;    % 


box_mp = [0.05*min(proteinMass):5e-15:1.25*max(proteinMass)];
box_mc = [0.5*min(cellMass):5e-15:1.25*max(cellMass)];
box_mu = [0*min(mu):0.001:1.1*max(mu)];
mp=[];
for k=1:length(box_mp)
  mp(k,1) = box_mp(k)/(c0+c1*box_mp(k));  
end
mc=[];
for k=1:length(box_mc)
  mc(k,1) = box_mc(k)/(d0+d1*box_mc(k));  
end

Fontsize=16;

f1=figure(1)
subplot(121)   
plot(box_mp,mp,'b-','LineWidth',3),ylabel('\mu'), xlabel('Protein DW')
hold on
plot(proteinMass,mu,'ks','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated \mu=m_p/(c_0+c_1*m_p)', 'Actual \mu','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off
subplot(122)   
plot(box_mc,mc,'b-','LineWidth',3),ylabel('\mu'), xlabel('Cell DW')
hold on
plot(cellMass,mu,'ks','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated \mu=m_c/(d_0+d_1*m_c)', 'Actual \mu','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off

f2=figure(2);
subplot(121)
plot(box_mu, box_mu./(alpha+beta_m*box_mu),'b-','LineWidth',3)
hold on
plot(box_mu, box_mu./(alpha_Z+beta_Z_m*box_mu),'g-','LineWidth',3)
plot(mu,1./(C_exp + D_exp),'sk','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Our estimated \mu/(\alpha+\beta*\mu)','Zheng"s estimated \mu/(\alpha+\beta*\mu)', '1/(C+D)','Location','southeast')
xlabel('Specific growth rate, \mu (min^{-1})')
ylabel('\mu/(\alpha+\beta*\mu), 1/(C+D)')
ax = gca;
ax.FontSize = 16;
hold off
subplot(122)
plot(box_mc,mc,'b-','LineWidth',3),ylabel('Growth rate \mu (min^{-1})'), xlabel('Cell DW (g)')
hold on
m0=545e-15; % value from Zheng's paper
%m0=112e-15; % our estimation in Main_fit_Fengwei_Zheng_mu.m
plot(box_mc,(box_mc-m0*alpha)/m0/beta_m,'g-','LineWidth',3 )
plot(box_mc,(box_mc-m0*alpha_Z)/m0/beta_Z_m,'r-','LineWidth',3 )

plot(cellMass,mu,'ks','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated \mu=m_c/(d_0+d_1*m_c)','Estimated \mu=(m_c-\alpha m_0)/m_0/\beta','Estimated \mu=(m_c-\alpha_Z m_0)/m_0/\beta_Z', 'Actual \mu','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off

% %f1.WindowState = 'maximized';
% %exportgraphics(f1,'../Texte/images/f_CpD_mu.png','Resolution',300)

f3=figure(3);
subplot(121)
plot(box_mc,mc,'b-','LineWidth',3),ylabel('Growth rate \mu (min^{-1})'), xlabel('Cell DW (g)')
hold on
m0=545e-15; % value from Zheng's paper
%m0=112e-15; % our estimation in Main_fit_Fengwei_Zheng_mu.m
%plot(box_mc,(box_mc-m0*alpha)/m0/beta_m,'g-','LineWidth',3 )
%plot(box_mc,(box_mc-m0*alpha_Z)/m0/beta_Z_m,'c-','LineWidth',3 )
plot(box_mc,log(box_mc/123.022e-15)/68.5547,'r-','LineWidth',3 )
plot(box_mc,(box_mc+1.08179e-13)/3.28694e-11,'k-','LineWidth',3 )
plot(cellMass,mu,'ks','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('\mu=m_c/(d_0+d_1*m_c)','Zheng best fit','Zheng original fit','Fengwei-Zheng best fit','Linear fit',  'Actual \mu','Location','southeast')
legend('\mu=m_c/(d_0+d_1*m_c)','Fengwei-Zheng best fit','Linear fit',  'Actual \mu','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off
subplot(122)
plot(mu,cellMass,'ks','MarkerSize',15,'MarkerFaceColor','k'), ylabel('m_c (DW)'), xlabel('Growth rate \mu (min^{-1})')
hold on
grid on
plot(box_mu,123.022e-15*exp(68.5547*box_mu),'b-','LineWidth',3)
plot(box_mu,1.9535e-16*exp(15.2264*box_mu.^0.16788),'c-','LineWidth',3)
plot(box_mu,-1.08179e-13 + 3.28694e-11*box_mu,'k-','LineWidth',3)

legend('Measured cell DW', 'm0*exp(\gamma*\mu)', 'm0*exp(\gamma*\mu^{\rho})', 'a+b*\mu','Location','southeast')
ax = gca;
ax.FontSize = 16;
hold off


f4=figure(4);
subplot(121)
plot(box_mu, box_mu./(alpha+beta_m*box_mu),'b-','LineWidth',3)
hold on
%plot(box_mu, box_mu./(alpha_Z+beta_Z_m*box_mu),'g-','LineWidth',3)
plot(mu,1./(C_exp + D_exp),'sk','MarkerSize',15,'MarkerFaceColor','k')
grid on
legend('Estimated $\mu/(\alpha+\beta\mu)$', 'Data $1/(C+D)$','Location','southeast','FontSize',Fontsize,'Interpreter','latex')
xlabel('Specific growth rate, $\mu (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
ylabel('$\mu/(\alpha+\beta\mu)$, $1/(C+D)$','FontSize',Fontsize,'Interpreter','latex')
ax = gca;
ax.FontSize = 16;
hold off
subplot(122)
plot(alpha+beta_m*mu,cellMass*1e15,'ks','MarkerSize',15,'MarkerFaceColor','k'), ylabel('$m_{cDW}\, (fg)$','FontSize',Fontsize,'Interpreter','latex')
xlabel('$\mu (C+D)$','FontSize',Fontsize,'Interpreter','latex')
hold on
grid on
plot(alpha+beta_m*box_mu,1e15*123.022e-15*exp(68.5547*box_mu),'b-','LineWidth',3)
plot(alpha+beta_m*box_mu,1e15*1.9535e-16*exp(15.2264*box_mu.^0.16788),'c-','LineWidth',3)
plot(alpha+beta_m*box_mu,1e15*(-1.08179e-13 + 3.28694e-11*box_mu),'k-','LineWidth',3)

legend('Measured cell DW (fg)', '$m_0e^{\gamma\mu}$', '$m_0 e^{\gamma\mu^{\rho}}$', '$a+b \mu$','Location','southeast','FontSize',Fontsize,'Interpreter','latex')
ax = gca;
ax.FontSize = 16;
hold off

%f4.WindowState = 'maximized';
% %exportgraphics(f4,'../images_def/﻿f_mc_muCpD.png','Resolution',300)


% % Our best estimated parameters using MATLAB script main_fit_CpD_mu.m
% % for mu/(alpha+beta*mu)=1/(C+D)
% 	alpha= 0.277079;
% 	beta_h= 0.912718; %(hours)
%     beta_m= beta_h*60; %(mins)
%     
%     % Our best estimated parameters using MATLAB script main_fit_Fengwei_Zheng_mu.m
% % for mc = m0*exp(gamma*mu):
%     m0 = 123.022e-15; %(g)
%     gamma = 68.5547; (mins)   % 

%%% RESULTS OBTAINED MODEL THREE PARAMETERS:mc = m0*exp(gamma*mu^rho)
% Best solution value		2.50967e-29
% Decision vector
% 	m0= 1.9535e-16
% 	gamma= 15.2264
% 	rho= 0.16788

%%% RESULTS OBTAINED MODEL TWO PARAMETERS: mc = a+b*mu
% Best solution value		1.23451e-27
% Decision vector
% 	-1.08179e-13
% 	3.28694e-11