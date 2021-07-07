
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               
%
% Evaluation of the function YpmRNA = b*x/(1+\mu_m*x)
% 
% with b= 0.62*nu/le,       x= f(si)*KC0*r/dm
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

%%%%%%% GENERAL CELL DATA:

mu_Hausser= log(2)/21.5; %0.032
td=20:5:60;
mu = log(2)./td;
le=[18,25,31];
nu_eff=0.62*1260./le;
general_dm= log(2)/2.5;
long_dm= log(2)/7.5;
 
Fontsize=16;

f1=figure(1)
subplot(121)
fplot(@(x)nu_eff(1).*x./general_dm./(1+mu_Hausser./general_dm.*x),[0 300],'Linewidth',3)
hold on
fplot(@(x)nu_eff(2).*x./general_dm./(1+mu_Hausser./general_dm.*x),[0 300],'Linewidth',3)
fplot(@(x)nu_eff(3).*x./general_dm./(1+mu_Hausser./general_dm.*x),[0 300],'Linewidth',3)
%fplot(@(x)nu_eff.*x,[0 1000],'Linewidth',3,'Color','r')
%fplot(@(x)nu_eff./mu_Hausser,[0 1000],'Linewidth',3,'Color','g')
grid on
xlabel('$K_C^0 r$','FontSize',Fontsize,'Interpreter','latex'), ylabel('$Y_{p/mRNA},\,\,\, d_m=\log(2)/2.5$','FontSize',Fontsize+1,'Interpreter','latex')
legend('$l_e=18$','$l_e=25$','$l_e=31$','FontSize',Fontsize,'Interpreter','latex', 'Location','southeast')

ax = gca;
 ax.YLim = [0 1400];
ax.FontSize = 13;
subplot(122)
fplot(@(x)nu_eff(1).*x./long_dm./(1+mu_Hausser./long_dm.*x),[0 300],'Linewidth',3)
hold on
fplot(@(x)nu_eff(2).*x./long_dm./(1+mu_Hausser./long_dm.*x),[0 300],'Linewidth',3)
fplot(@(x)nu_eff(3).*x./long_dm./(1+mu_Hausser./long_dm.*x),[0 300],'Linewidth',3)
%fplot(@(x)nu_eff.*x,[0 1000],'Linewidth',3,'Color','r')
%fplot(@(x)nu_eff./mu_Hausser,[0 1000],'Linewidth',3,'Color','g')
grid on
xlabel('$K_C^0 r$','FontSize',Fontsize,'Interpreter','latex'), ylabel('$Y_{p/mRNA},\,\,\, d_m=\log(2)/7.5$','FontSize',Fontsize+1,'Interpreter','latex')
legend('$l_e=18$','$l_e=25$','$l_e=31$','FontSize',Fontsize,'Interpreter','latex', 'Location','southeast')
ax = gca;
 ax.YLim = [0 1400];
ax.FontSize = 13;
%print -dpng ../Texte/images/f_YpmRNA_RBS.png
%f1.WindowState = 'maximized';
%exportgraphics(f1,'../images_def/f_YpmRNA_RBS.png','Resolution',300)
 hold off

figure(2)
fplot(@(x)nu_eff(2).*x./general_dm./(1+mu_Hausser./general_dm.*x),[0 60],'Linewidth',3)
hold on 
fplot(@(x)nu_eff(2).*x./general_dm,[0 60],'Linewidth',3,'Color','g')
grid on
xlabel('K_C^0 r'), ylabel('YpmRNA, d_m=log(2)/2.5')
ax = gca;
 ax.YLim = [0 1000];
ax.FontSize = 13;

general_dm/mu_Hausser/0.2

long_dm/mu_Hausser/0.2


