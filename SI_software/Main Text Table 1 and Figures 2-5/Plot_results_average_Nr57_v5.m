%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp5.m for Nr=57
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

Bremer_exp_data.Phi_t = 0.7796; %'Fraction of mature available ribosomes relative to their total number   (adim)   

%%%%%%% General values of cell parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

model_p.l_e = 25; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.nu_max = 1260; % Maximum effective translation rate per ribosome (aa/min)
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57; %Number of protein types explaining 99% J's in Hausser data (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)		
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa)
Da = 1.6605e-24; % g/Da
model_p.ribosome_mass = 2.29*model_p.m_aa*model_p.N_r*model_p.lp_r; % mp_estimated includes the RNA component of ribosomes.
%model_p.model_mass=1; % Protein mass model  m =m0exp(βμ)    m0=77.3748e−15(g)    β=61.7813(min)
model_p.model_mass=2; % Protein mass model  m =m0exp(βμ^gamma)  

%%%%% Weighted average values of estimated parameters

model_p.ku_r =  130.7432; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = 3.1120; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  5.8138; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  11.8005; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = 5.3606; % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = 0.0278; % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.Phi_t = 0.7774; % Average ratio between mature available and total number of ribosomes 

best_N57_params=[model_p.ku_r,model_p.ku_nr, model_p.kb_r, model_p.kb_nr, model_p.Omega_r, model_p.Omega_nr,model_p.Phi_t ];

model_p.mu_estimated = Bremer_exp_data.mu; %Initial estimations for the function eval_Mu_ss_wildtype_mp5.m
model_p.Phi_b = ones(5,1); 
model_p.Phi_r = ones(5,1); 

%%%%% EVALUATION OF THE AVERAGE ESTIMATION:

evaluate_results=true;
if evaluate_results==true

[ke,KC0_r,KC0_nr, mu_estimated_profile,mp_estimated,mu_r_profile] = eval_Mu_ss_wildtype_mp5(best_N57_params); 

JSum=model_p.Phi_b./(1-model_p.Phi_b);
JNr=model_p.Phi_r.*(1+JSum);
Fraction_r=JNr./JSum;
varphi_R = Bremer_exp_data.rt*model_p.ribosome_mass./Bremer_exp_data.mp;
varphi_P = 1-varphi_R;
estimated_varphi_R = Fraction_r;
estimated_varphi_P = 1-estimated_varphi_R;

r_estimated = mu_r_profile./mu_estimated_profile;

rt_estimated = (1+JSum).*r_estimated/model_p.Phi_t;

BetaP_R_estimated = JNr.* r_estimated.*Bremer_exp_data.nu_t .*model_p.dm_r./model_p.lp_r./model_p.Omega_r/model_p.N_r % Translation rate beta_p
BetaP_P_estimated = (JSum-JNr).* r_estimated.*Bremer_exp_data.nu_t.*model_p.dm_nr./model_p.lp_nr./model_p.Omega_nr/model_p.N_nr
Limit_r=14*model_p.Omega_r
Limit_nr=14*model_p.Omega_nr

figures_plot=true;

if figures_plot==true
 
f1= figure(1)
subplot(131)
plot(Bremer_exp_data.mu, mu_estimated_profile,'b*','MarkerSize',10,'LineWidth',3), xlabel('estimated \mu (min^{-1})'), ylabel('experimental \mu (min^{-1})'), grid on
hold on
fplot(@(x)x, [min(Bremer_exp_data.mu) max(Bremer_exp_data.mu)],'r-','Linewidth',2)
ax = gca;
ax.FontSize = 13;
hold off
subplot(132)
plot(mu_estimated_profile, r_estimated,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('r (molec)'), grid on
hold on
%plot(Bremer_exp_data.mu, 350*Bremer_exp_data.rt/Bremer_exp_data.rt(5),'r*','MarkerSize',10,'LineWidth',3)
%legend('Estimated free ribosomes', 'Initial estimation','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
subplot(133)
plot(mu_estimated_profile, rt_estimated,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('r_t (molec)'), grid on
hold on
plot(Bremer_exp_data.mu, Bremer_exp_data.rt ,'r*','MarkerSize',10,'LineWidth',3)
legend('Estimated r_t', 'Experimental  r_t','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off

%3


f2= figure(2)
subplot(131)
plot(mu_estimated_profile, ke,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('k_e'), grid on
ax = gca;
ax.FontSize = 13;
hold off
subplot(132)
plot(mu_estimated_profile,KC0_r,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('K_{C0}^r, K_{C0}^{nr}'), grid on
hold on
plot(mu_estimated_profile, KC0_nr,'r*','MarkerSize',10,'LineWidth',3)
legend('K_{C0}^r', 'K_{C0}^{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
subplot(133)
plot(log(mu_r_profile),model_p.dm_r./KC0_r,'b*','MarkerSize',10,'LineWidth',3), xlabel('log(\mu r) (free ribosomes\cdot min^{-1})'), ylabel('d_{m,k}/K_{C0}^k'), grid on
hold on
plot(log(mu_r_profile), model_p.dm_nr./KC0_nr,'r*','MarkerSize',10,'LineWidth',3)
legend('d_{m,r}/K_{C0}^r', 'd_{m,nr}/K_{C0}^{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
%exportgraphics(f2,'./RBS_fractions_v5.png','Resolution',300)

 f3=figure(3);
plot(Bremer_exp_data.mu, varphi_R,'bd','MarkerSize',10,'LineWidth',3), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'rd','MarkerSize',10,'LineWidth',3)
plot(mu_estimated_profile, estimated_varphi_R,'bo','MarkerSize',10,'LineWidth',3) 
plot(mu_estimated_profile, estimated_varphi_P,'ro','MarkerSize',10,'LineWidth',3)

legend('Exp. \phi_R','Exp. \phi_P','Extimated \phi_R','Estimated. \phi_P')
 xlabel('\mu (min^{-1})'), ylabel('\phi_R, \phi_P')
 title(' Ribosomal (\phi_R) and non-ribosomal protein  (\phi_P) mass fractions')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off
%exportgraphics(f3,'./mass_fractions_v5.png','Resolution',300)

global data_log10mur;
global data_log10r;
data_log10mur=log(mu_r_profile);
data_log10r=log(r_estimated);
problem2.f='mse_wildtype';    %script with the cost function to optimize
problem2.x_L=[-5,0]; % minimum expected values for ku_r, ku_nr, kb_r, kb_nr,omega_r, omega_nr, Phi_t
problem2.x_U=[5,1.5]; %  maximum expected values, Phi_t
opts2.maxeval=7500; 
Results2=MEIGO(problem2,opts2,'ESS'); 
% Best solution value		0.0752833
% Decision vector
% 	4.10928
% 	0.792795
 f4=figure(4);
 R_s=[model_p.Phi_t*rt_estimated' ;r_estimated'];
 subplot(121)
area(mu_estimated_profile,R_s','LineWidth',3), grid on
legend('Mature available ribosomes','Free ribosomes')
 xlabel('\mu (min^{-1})'), ylabel('number ribosomes')
 title(' Estimated ribosomes')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off
  subplot(122)
  plot(log(mu_r_profile),log(r_estimated),'b*','MarkerSize',10,'LineWidth',3)
  hold on
  fplot(@(x)Results2.xbest(1,1) + Results2.xbest(1,2)*x,[-4,4],'r-','LineWidth',3)
  xlabel('log(\mu r) (free ribosomes\cdot min^{-1})'), ylabel('log(r)')
  grid on
 ax = gca;
 ax.FontSize = 13;
hold off
%exportgraphics(f4,'./r_ra_sums_v5.png','Resolution',300)


f6= figure(6)
subplot(131)
plot(mu_estimated_profile, JNr,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('N_xJ_x'), grid on
hold on
plot(mu_estimated_profile, JSum-JNr,'r*','MarkerSize',10,'LineWidth',3)
legend('N_rJ_r', 'N_{nr}J_{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
subplot(132)
plot(mu_r_profile,JNr,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu r (free ribosomes\cdot min^{-1}'), ylabel('N_xJ_x'), grid on
hold on
plot(mu_r_profile,JSum-JNr,'r*','MarkerSize',10,'LineWidth',3)
legend('N_rJ_r', 'N_{nr}J_{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
subplot(133)
plot(log(mu_r_profile),JNr,'b*','MarkerSize',10,'LineWidth',3),
xlabel('log(\mu r) (free ribosomes\cdot min^{-1})'), ylabel('N_xJ_x'), grid on
hold on
plot(log(mu_r_profile), JSum-JNr,'r*','MarkerSize',10,'LineWidth',3)
legend('N_rJ_r', 'N_{nr}J_{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off
%exportgraphics(f6,'./﻿J_fractions_mur_v5c.png','Resolution',300)

f7=figure(7);

subplot(121)
plot(Bremer_exp_data.mu, varphi_R,'bd','MarkerSize',10,'LineWidth',3), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'rd','MarkerSize',10,'LineWidth',3)
plot(mu_estimated_profile, estimated_varphi_R,'bo','MarkerSize',10,'LineWidth',3) 
plot(mu_estimated_profile, estimated_varphi_P,'ro','MarkerSize',10,'LineWidth',3)

legend('Exp. \phi_R','Exp. \phi_P','Extimated \phi_R','Estimated. \phi_P')
 xlabel('\mu (min^{-1})'), ylabel('\phi_R, \phi_P')
 title(' Ribosomal (\phi_R) and non-ribosomal protein  (\phi_P) mass fractions')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off

subplot(122)
plot(Bremer_exp_data.mu, Bremer_exp_data.mp.*varphi_R,'bd','MarkerSize',10,'LineWidth',3), grid on
hold on
plot(Bremer_exp_data.mu, Bremer_exp_data.mp.*varphi_P,'rd','MarkerSize',10,'LineWidth',3)
plot(Bremer_exp_data.mu, Bremer_exp_data.mp,'kd','MarkerSize',10,'LineWidth',3)
plot(mu_estimated_profile, mp_estimated.*estimated_varphi_R,'bo','MarkerSize',10,'LineWidth',3) 
plot(mu_estimated_profile, mp_estimated.*estimated_varphi_P,'ro','MarkerSize',10,'LineWidth',3)
plot(mu_estimated_profile,mp_estimated,'ko','MarkerSize',10,'LineWidth',3)

legend('Exp. m_{R}','Exp. m_{P}','Exp. m_{R+P}','Extimated m_R','Estimated m_P','Estimated m_{R+P}')
 xlabel('\mu (min^{-1})'), ylabel('m_x')
 title('Total, ribosomal and non-ribosomal protein masses')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off


end
end



