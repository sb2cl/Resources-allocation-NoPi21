%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp_v6.m for Nr=57
%
%  This script runs the average model obtained from the best 25 runs and
%  plots the results.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mex;
clear all;
close all;
warning off;
global model_p;
global Bremer_exp_data;

%  We use the values in Bremer (cite) to compare with
Bremer_exp_data.mu = [
    log(2)/100
    log(2)/60
    log(2)/40
    log(2)/30
    log(2)/24]; %1/min

Bremer_exp_data.mp = [
    100
    156
    234
    340
    450]*1e-15; % grams

% Effective translation rate. aa/min
Bremer_exp_data.nu_t = [
    12
    16
    18
    20
    21]*60;

% Number of ribosomes (mature and inmature). molecs
Bremer_exp_data.rt = [
    6800
    13500
    26300
    45100
    72000];

Bremer_exp_data.Phi_t = 0.7796; %'Fraction of active ribosomes relative to their total number   (adim)   

%%%%%%% General values of cell parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

model_p.l_e = 25; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.nu_max = 1260; % Maximum effective translation rate per ribosome (aa/min)
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57; %Number of protein types explaining 99% J's in Hausser data (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)		
model_p.WEm_r = 1.2891; % Average weight (1+1/Emr) for ribosomal protein-coding genes
model_p.WEm_nr =  1.1575; % Average weight (1+1/Emp) for non-ribosomal protein-coding genes
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa)
Da = 1.6605e-24; % g/Da
model_p.ribosome_mass = 2.29*model_p.m_aa*model_p.N_r*model_p.lp_r; % mp_estimated includes the RNA component of ribosomes.
%model_p.model_mass=1; % Protein mass model  m =m0exp(βμ)    m0=77.3748e−15(g)    β=61.7813(min)
model_p.model_mass=2; % Protein mass model  m =m0exp(βμ^gamma)  

%%%%% Weighted average values of estimated parameters

model_p.ku_r =  133.04; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = 3.09; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  5.88; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  13.41; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = 5.72; % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = 0.028; % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.Phi_t = 0.969; % Average ratio between mature available and total number of ribosomes 

best_N57_params=[model_p.ku_r,model_p.ku_nr, model_p.kb_r, model_p.kb_nr, model_p.Omega_r, model_p.Omega_nr,model_p.Phi_t ];

model_p.mu_estimated = Bremer_exp_data.mu; %Initial estimations to be used by the function eval_Mu_ss_wildtype_mp_v6.m
model_p.Phi_b = 0.99*ones(5,1); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_b_t = 0.83*ones(5,1); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_r_t = 0.3*ones(5,1);  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_p_t = 0.55*ones(5,1);  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)

%%%%% EVALUATION OF THE AVERAGE ESTIMATION:

evaluate_results=true;
if evaluate_results==true
    
    Weighted_average_results_Nr57=[model_p.ku_r;model_p.ku_nr;model_p.kb_r;model_p.kb_nr;model_p.Omega_r;model_p.Omega_nr;model_p.Phi_t];
    [ke,KC0_r,KC0_nr, mu_estimated_profile,mp_estimated,mu_r_profile] = eval_Mu_ss_wildtype_mp_v6(Weighted_average_results_Nr57'); 

JWSum=model_p.Phi_b./(1-model_p.Phi_b); %Total weighted sum of J's
JNr=model_p.Phi_r_t.*(1+JWSum); %Sum of ribosomal J's
JNnr=model_p.Phi_p_t.*(1+JWSum); %Sum of non-ribosomal J's
JSum = JNr + JNnr; %Total sum of J's
Fraction_r=JNr./JSum; % Estimated fraction of ribosomal mass
varphi_R = Bremer_exp_data.rt*model_p.ribosome_mass./Bremer_exp_data.mp; %Experimental fraction of ribosomal mass
varphi_P = 1-varphi_R;
estimated_varphi_R = Fraction_r;
estimated_varphi_P = 1-estimated_varphi_R;

r_estimated = mu_r_profile./mu_estimated_profile; %r_estimated = estimated_flux_mur/estimated_mu
rt_estimated = (1+JWSum).*r_estimated/model_p.Phi_t; %Estimated total number of ribosomes: r_a = Phi_t*r_t -> r_t=r_a/Phi_t = r/(1+JWSum)/Phi_t
r_mature_available = (1+JWSum).*r_estimated;
r_active_estimated = (model_p.Phi_r_t+   model_p.Phi_p_t).*r_mature_available;  
r_bound_at_RBS=r_mature_available-r_active_estimated;
r_immature= rt_estimated-(r_mature_available+ r_estimated);

%Next we can estimate limits for translation rate as those encountered by
%Hausser:
%BetaP_R_estimated = JNr.* r_estimated.*Bremer_exp_data.nu_t .*model_p.dm_r./model_p.lp_r./model_p.Omega_r/model_p.N_r % Translation rate beta_p
%BetaP_P_estimated = (JSum-JNr).* r_estimated.*Bremer_exp_data.nu_t.*model_p.dm_nr./model_p.lp_nr./model_p.Omega_nr/model_p.N_nr
%Limit_r=14*model_p.Omega_r
%Limit_nr=14*model_p.Omega_nr

figures_plot=true;

if figures_plot==true
 
f1= figure(1)
subplot(131)
plot(Bremer_exp_data.mu, mu_estimated_profile,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineStyle','none'), 
xlabel('estimated $\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('experimental  $\mu (\mathrm{min}^{-1})$','Interpreter','latex'), grid on
hold on
fplot(@(x)x, [min(Bremer_exp_data.mu) max(Bremer_exp_data.mu)],'r-','Linewidth',2,'LineStyle','-')
ax = gca;
ax.FontSize = 14;
hold off
subplot(132)
plot(mu_estimated_profile, r_estimated,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$r$ (molec)','Interpreter','latex'), grid on
hold on
%plot(Bremer_exp_data.mu, 350*Bremer_exp_data.rt/Bremer_exp_data.rt(5),'r*','MarkerSize',10,'LineWidth',3)
%legend('Estimated free ribosomes', 'Initial estimation','Location','southeast')
ax = gca;
ax.FontSize = 14;
hold off
subplot(133)
plot(mu_estimated_profile, rt_estimated,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$r_t$ (molec)','Interpreter','latex'), grid on
hold on
plot(Bremer_exp_data.mu, Bremer_exp_data.rt ,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--')
legend('Estimated $r_t$', 'Experimental  $r_t$','Location','southeast','FontSize',14,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
hold off
%exportgraphics(f1,'./estimation_fractions_v6.png','Resolution',300)


f2= figure(2)
subplot(131)
plot(mu_estimated_profile, ke,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$k_e\, (\mathrm{min}^{-1})$','Interpreter','latex'), grid on
ax = gca;
ax.FontSize = 14;
hold off
subplot(132)
plot(mu_estimated_profile,KC0_r,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$K_{C0}^r$, $K_{C^0}^{nr}\, (\mathrm{molec.}^{-1})$','Interpreter','latex'), grid on
hold on
plot(mu_estimated_profile, KC0_nr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--')
legend('$K_{C0}^r$', '$K_{C^0}^{nr}$','Interpreter','latex','Location','northeast','FontSize',16)
ax = gca;
ax.FontSize = 14;
hold off
subplot(133)
plot(log(mu_r_profile),model_p.dm_r./KC0_r,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\log{\mu r}\, (\mathrm{molec.} \cdot \mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$d_{m,k}/K_{C^0}^k\, (\mathrm{molec.} \cdot \mathrm{min}^{-1})$','Interpreter','latex'), grid on
hold on
plot(log(mu_r_profile), model_p.dm_nr./KC0_nr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--')
legend('$d_{m,r}/K_{C^0}^r$', '$d_{m,nr}/K_{C^0}^{nr}$','Location','southeast','Interpreter','latex','FontSize',16)
ax = gca;
ax.FontSize = 14;
hold off
%exportgraphics(f2,'./RBS_fractions_v6.png','Resolution',300)

 f3=figure(3);
plot(Bremer_exp_data.mu, varphi_R,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
plot(mu_estimated_profile, estimated_varphi_R,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b') 
plot(mu_estimated_profile, estimated_varphi_P,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')

legend('Measured $ \phi_R$','Measured $\phi_P$','Extimated $\phi_R$','Estimated. $\phi_P$','Interpreter','latex','Location','southeast','FontSize',16)
 xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$\phi_R$, $\phi_P$','Interpreter','latex')
 title(' Ribosomal ($\phi_R$) and non-ribosomal  ($\phi_P$) protein mass fractions','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off

global data_log10mur;
global data_log10r;
data_log10mur=log(mu_r_profile);
data_log10r=log(r_estimated);
problem2.f='mse_wildtype';    %script with the cost function to optimize
problem2.x_L=[-5,0]; % minimum expected values for ku_r, ku_nr, kb_r, kb_nr,omega_r, omega_nr, Phi_t
problem2.x_U=[5,1.5]; %  maximum expected values, Phi_t
opts2.maxeval=7500; 
Results2=MEIGO(problem2,opts2,'ESS'); 
% Best solution value		0.0603
% Decision vector
%   4.069 
% 	0.782


 f4=figure(4);
 subplot(131)
  R_s=[r_active_estimated'; r_bound_at_RBS' ;r_estimated';r_immature'];
area(mu_estimated_profile,R_s','LineWidth',2), grid on
legend('Active bound ribosomes', 'RBS bound ribosomes','Free ribosomes','Immature ribosomes', 'Location','northwest','FontSize',16)
 xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('number of ribosomes')
 title(' Estimated ribosomes','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off
 subplot(132)
 R_s=[r_active_estimated'./rt_estimated'; r_bound_at_RBS'./rt_estimated' ;r_estimated'./rt_estimated';r_immature'./rt_estimated'];
area(mu_estimated_profile,R_s','LineWidth',2), grid on
legend('Active bound ribosomes', 'RBS bound ribosomes','Free ribosomes','Immature ribosomes', 'Location','southeast','FontSize',16)
 xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('fraction of ribosomes')
 title(' Estimated fractions of ribosomes','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off 
  subplot(133)
  plot(log(mu_r_profile),log(r_estimated),'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b')
  hold on
  fplot(@(x)Results2.xbest(1,1) + Results2.xbest(1,2)*x,[-4,4],'LineWidth',2,'LineStyle','-','Color','r')
  xlabel('$\log{mu r}$ (ribosomes$\cdot \mathrm{min}^{-1}$)','Interpreter','latex'), ylabel('$\log{r}$','Interpreter','latex')
  grid on
 ax = gca;
 ax.FontSize = 14;
hold off
%exportgraphics(f4,'../images_def/r_ra_sums_v6.png','Resolution',300)


f6= figure(6)
% subplot(131)
subplot(121)
plot(mu_estimated_profile, JNr,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'), 
xlabel('$\mu$ $(\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$N_xJ_x$','Interpreter','latex'), grid on
hold on
plot(mu_estimated_profile, JNnr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
legend('$N_rJ_r$', '$N_{nr}J_{nr}$','Location','southeast','Interpreter','latex','FontSize',16)
ax = gca;
ax.FontSize = 14;
hold off
% subplot(132)
% plot(mu_r_profile,JNr,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'), 
% xlabel('$\mu r$ (ribosomes$\cdot \mathrm{min}^{-1}$)','Interpreter','latex'), ylabel('$N_xJ_x$','Interpreter','latex'), grid on
% hold on
% plot(mu_r_profile,JNnr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
% legend('$N_rJ_r$', '$N_{nr}J_{nr}$','Location','southeast','Interpreter','latex','FontSize',16)
% ax = gca;
% ax.FontSize = 14;
% hold off
% subplot(133)
subplot(122)
plot(log(mu_r_profile),JNr,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'),
xlabel('$\mu r$ (ribosomes$\cdot \mathrm{min}^{-1}$)','Interpreter','latex'), ylabel('$N_xJ_x$','Interpreter','latex'), grid on
hold on
plot(log(mu_r_profile), JNnr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
legend('$N_rJ_r$', '$N_{nr}J_{nr}$','Location','southeast','Interpreter','latex','FontSize',16)
ax = gca;
ax.FontSize = 14;
hold off
%exportgraphics(f6,'../﻿images_def/J_fractions_mur_v6.png','Resolution',300)


f7=figure(7);
subplot(121)
plot(Bremer_exp_data.mu, varphi_R,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
plot(mu_estimated_profile, estimated_varphi_R,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b') 
plot(mu_estimated_profile, estimated_varphi_P,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')

legend('Measured $ \phi_R$','Measured $\phi_P$','Extimated $\phi_R$','Estimated. $\phi_P$','Interpreter','latex','Location','southeast','FontSize',16)
 xlabel('$\mu\, (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$\phi_R$, $\phi_P$','Interpreter','latex')
 title(' Ribosomal ($\phi_R$) and non-ribosomal  ($\phi_P$) protein mass fractions','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off

subplot(122)
plot(Bremer_exp_data.mu, Bremer_exp_data.mp.*varphi_R,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b'), grid on
hold on
plot(Bremer_exp_data.mu, Bremer_exp_data.mp.*varphi_P,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
plot(Bremer_exp_data.mu, Bremer_exp_data.mp,'Marker','s','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','k')
plot(mu_estimated_profile, mp_estimated.*estimated_varphi_R,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b') 
plot(mu_estimated_profile, mp_estimated.*estimated_varphi_P,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r')
plot(mu_estimated_profile,mp_estimated,'Marker','o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','k')

legend('Measured $m_{R}$','Measured $m_{P}$','Measured $m_{R+P}$','Extimated $m_R$','Estimated $m_P$','Estimated $m_{R+P}$','Interpreter','latex','Location','southeast','FontSize',16)
 xlabel('$\mu\,  (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$m_x$','Interpreter','latex')
 title('Protein mass','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off
%exportgraphics(f7,'../images_def/mass_fractions_v6.png','Resolution',300)


end
end



