%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp5.m for Nr=57
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mex;
clear all;
close all;
warning off;
global model_p;
global Bremer_exp_data;

num_runs=200;

Best_Results_Nr57=[];
for k=1:num_runs
 main_fitMu_ss_wildtype_mp5
 Best_Results_Nr57= [Best_Results_Nr57, [Results.fbest, Results.xbest]'];
end
 
[const_index_sorted,cost_indexl_sorted_index]=sort(Best_Results_Nr57(1,:),'descend');


Best_Results_Nr57_sorted = Best_Results_Nr57(:,cost_indexl_sorted_index);

Best_Results_Nr57_sorted(1,:)= 1./Best_Results_Nr57_sorted(1,:);
% Best result is the last one.
% Get rid of  worst results. 
num_bad=175;
Best_Results_Nr57_sorted_choice=Best_Results_Nr57_sorted(:,num_bad+1:end);

Weighted_sum_results_Nr57 = zeros(size(Best_Results_Nr57_sorted_choice,1),1);
Sum_weighs_results_Nr57 = 0;

for k=1:size(Best_Results_Nr57_sorted_choice,2)
    Weighted_sum_results_Nr57 = Weighted_sum_results_Nr57 + Best_Results_Nr57_sorted_choice(:,k)*Best_Results_Nr57_sorted_choice(1,k);
    Sum_weighs_results_Nr57 = Sum_weighs_results_Nr57 + Best_Results_Nr57_sorted_choice(1,k);
end
Weighted_average_results_Nr57= Weighted_sum_results_Nr57/Sum_weighs_results_Nr57

% NOW we simulate the steady %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%% Weighted average values of estimated parameters

model_p.ku_r =  Weighted_average_results_Nr57(2,1); %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = Weighted_average_results_Nr57(3,1); %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  Weighted_average_results_Nr57(4,1); %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  Weighted_average_results_Nr57(5,1); %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = Weighted_average_results_Nr57(6,1); % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr =Weighted_average_results_Nr57(7,1); % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.Phi_t = Weighted_average_results_Nr57(8,1); % Average ratio between mature available and total number of ribosomes 

%%%%% EVALUATION OF THE AVERAGE ESTIMATION:

evaluate_results=true;
if evaluate_results==true

[ke,KC0_r,KC0_nr, mu_estimated_profile,mp_estimated,mu_r_profile] = eval_Mu_ss_wildtype_mp5(Weighted_average_results_Nr57(2:8,1)'); 

JSum=model_p.Phi_b./(1-model_p.Phi_b);
JNr=model_p.Phi_r.*(1+JSum);
Fraction_r=JNr./JSum;
varphi_R = Bremer_exp_data.rt*model_p.ribosome_mass./Bremer_exp_data.mp;
varphi_P = 1-varphi_R;
estimated_varphi_R = Fraction_r;
estimated_varphi_P = 1-estimated_varphi_R;

r_estimated = mu_r_profile./mu_estimated_profile;

rt_estimated = (1+JSum).*r_estimated/model_p.Phi_t;

 

figures_plot=true;

if figures_plot==true
 
f1= figure(1)
subplot(131)
plot(Bremer_exp_data.mu, mu_estimated_profile,'b*','MarkerSize',10,'LineWidth',3), xlabel('estimated \mu (min^{-1})'), ylabel('exp \mu (min^{-1})'), grid on
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

%exportgraphics(f1,'./estimation_fractions.png','Resolution',300)


f2= figure(2)
subplot(121)
plot(mu_estimated_profile, ke,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('k_e'), grid on
ax = gca;
ax.FontSize = 13;
hold off
subplot(122)
plot(mu_estimated_profile,KC0_r,'b*','MarkerSize',10,'LineWidth',3), xlabel('\mu (min^{-1})'), ylabel('K_{C0}^r, K_{C0}^{nr}'), grid on
hold on
plot(mu_estimated_profile, KC0_nr,'r*','MarkerSize',10,'LineWidth',3)
legend('K_{C0}^r', 'K_{C0}^{nr}','Location','southeast')
ax = gca;
ax.FontSize = 13;
hold off

%exportgraphics(f2,'./RBS_fractions.png','Resolution',300)

 f3=figure(3);
plot(Bremer_exp_data.mu, varphi_R,'rd','MarkerSize',10,'LineWidth',3), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'bd','MarkerSize',10,'LineWidth',3)
plot(mu_estimated_profile, estimated_varphi_R,'ro','MarkerSize',10,'LineWidth',3) 
plot(mu_estimated_profile, estimated_varphi_P,'bo','MarkerSize',10,'LineWidth',3)

legend('Exp. \phi_R','Exp. \phi_P','Extimated \phi_R','Estimated. \phi_P')
 xlabel('\mu (min^{-1})'), ylabel('\phi_R, \phi_P')
 title(' Ribosomal (\phi_R) and non-ribosomal protein  (\phi_P) mass fractions')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off
%exportgraphics(f3,'./mass_fractions.png','Resolution',300)

 f4=figure(4);
 R_s=[model_p.Phi_t*rt_estimated' ;r_estimated'];
area(mu_estimated_profile,R_s','LineWidth',3), grid on
legend('Mature available ribosomes','Free ribosomes')
 xlabel('\mu (min^{-1})'), ylabel('number ribosomes')
 title(' Estimated ribosomes')
 grid on
 ax = gca;
 ax.FontSize = 13;
hold off
 
%exportgraphics(f4,'./r_ra_sums.png','Resolution',300)

end
end

% Weighted_average_results_Nr57 =  mean of best 25 out of 200
% 
%   345.1734
%   130.7432
%     3.1120
%     5.8138
%    11.8005
%     5.3606
%     0.0278
%     0.7774
