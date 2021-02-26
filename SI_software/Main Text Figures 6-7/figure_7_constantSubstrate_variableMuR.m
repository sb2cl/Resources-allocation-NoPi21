%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state interactions host-protein considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp5.m for Nr=57
%   as analysed in >Analysis_resultys_estim_Nr57_v5.m and Plot_results_average_Nr57_v5.m
%   for the host dynamics
%  Here we get the values of productivity as a function of promoter and RBS
%  strengths for the protein of interest (A) for varying substrate
%  availability. We assume very low J's, as corresponding to an endogenous
%  protein
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear mex;
clear all;
close all;
warning off;

global model_p;

%  SUBSTRATE: We use the values of f(s_i) that include Bremer's ones as reference
% Obtained using Main_fit_si_mu.m
%f(s_i) = s_i/(K_s+s_i)
td_points=24:6:100;
Bremer_exp_data.mu= log(2)./td_points;
b0=4.45687e-05;
b1=19.9843;
a1=34.9866;
gamma1= 0.679116;
gamma2=1.06627;
Bremer_exp_data.f_si= (b0 + b1*Bremer_exp_data.mu.^gamma1)./(1+a1*Bremer_exp_data.mu.^gamma2);

f_substrate = Bremer_exp_data.f_si;

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

%%%%% Weighted average values of estimated parameters for the host

model_p.ku_r =  130.7432; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = 3.1120; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)

model_p.kb_r =  5.8138; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  11.8005; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)

model_p.Omega_r = 5.3606; % Average transcription rate for ribosomal proteins  (1/min)
model_p.Omega_nr = 0.0278; % Average transcription rate for non-ribosomal proteins  (1/min)
model_p.Phi_t = 0.7774; % Average ratio between mature available and total number of ribosomes

model_p.mu_estimated = Bremer_exp_data.mu; %Initial estimations for the function eval_host_protein_ss.m
model_p.Phi_b = ones(1,length(model_p.mu_estimated));
model_p.Phi_r = ones(1,length(model_p.mu_estimated));

%%%%%%% Initial values for the promoter and RBS strengths of the exogenous
%%%%%%% protein A

model_p.ku_A = 3.1120; %Dissotiation rate RBS-ribosome  (1/min)          Assume non-ribosomal protein
model_p.kb_A =  11.8005; %Assotiation rate RBS-ribosome  (1/min/molec)
model_p.Omega_A = 0.03; % Transcription rate   (1/min)
model_p.lp_A = model_p.lp_nr; % Assume protein A with mean characteristics
model_p.dm_A = model_p.dm_nr;

proteinA_params =[model_p.ku_A, model_p.kb_A, model_p.Omega_A];

%%%%% EVALUATION FOR AVERAGE RBS AND VARYING PROMOTER:

% varying_omega_A=[0.1,1,5,10,15]*model_p.Omega_A;
varying_omega_A=(0.0:0.1:1)*model_p.Omega_A*1;

%%%%% EVALUATION FOR AVERAGE  PROMOTER  AND VARYING RBS:
% Recall average values model_p.ku_A = 3.1120; %Dissotiation rate RBS-ribosome  (1/min)
% model_p.kb_A =  11.8005; %Assotiation rate RBS-ribosome  (1/min/molec)
range_Ku=[3,135];
range_Kb=[3,15];
range_Ku=[3,135];
range_Kb=[3,15];
n_levels=10;
discrete_step_Ku=round((max(range_Ku)-min(range_Ku))/n_levels,0);
discrete_step_Kb=round((max(range_Kb)-min(range_Kb))/n_levels,1);
dFF = fullfact([n_levels+1, n_levels+1]); % Full factorial experiment
varying_matrix_Ku_Kb_A = [3+(dFF(:,1)-1)*discrete_step_Ku, 3+(dFF(:,2)-1)*discrete_step_Kb];
varying_ratio_Kb_over_Ku_A= varying_matrix_Ku_Kb_A(:,2)./(varying_matrix_Ku_Kb_A(:,1)+model_p.nu_max/model_p.l_e); %RBS strength for substrate saturation
[sorted_Kb_over_Ku_A,index_sorted_Kb_over_Ku_A]=sort(varying_ratio_Kb_over_Ku_A,'ascend'); %sorted RBS strengths
sorted_varying_matrix_Ku_Kb_A = varying_matrix_Ku_Kb_A(index_sorted_Kb_over_Ku_A,:); % sorted matrix of (Ku,Kb) pairs according to ascending RBS strength for saturated substrate

%%%%%%% Not needed to check, for there is a monotonous relationship, but gives us the
%%%%%%% set of all RBS strengths for all values of the substrate concentration function
%%%%%%% f(s_i) needed for plotting
check_varying_ratio_Kb_over_Ku_A=[];
for k=1:length(f_substrate)
    check_varying_ratio_Kb_over_Ku_A(:,k)= sorted_varying_matrix_Ku_Kb_A(:,2) ./ ( sorted_varying_matrix_Ku_Kb_A(:,1) + model_p.nu_max/model_p.l_e*f_substrate(k) );
end
% check_varying_ratio_Kb_over_Ku_A
%%%%%%%%%%%%%%%%%%%%%%%%%

%varying_omega_A=[0.1,1,10,50,100,500,1000,1500,2500,5000,7500,10000]*model_p.Omega_A;
for p=1:length(varying_omega_A)
    for k=1:size(sorted_varying_matrix_Ku_Kb_A,1)
        % Recall average model_p.Omega_A = 0.03; % Transcription rate   (1/min)
        % Recall: varying_omega_A=[0.1, 1,10,50,100,500,1000,1500,3000,5000]*model_p.Omega_A;
        proteinA_params =[sorted_varying_matrix_Ku_Kb_A(k,1),sorted_varying_matrix_Ku_Kb_A(k,2),varying_omega_A(p)];
        
        [ke,KC0_r,KC0_nr, KC0_A,mu_estimated_profile,mp_estimated,mu_r_profile, JSum, JNr, J_A]  = eval_host_protein_ss(proteinA_params, f_substrate);
        
        estimated_fraction_A = J_A./JSum;
        estimated_fraction_R = JNr./JSum;
        estimated_fraction_P = 1-estimated_fraction_A-estimated_fraction_R;
        
        estimated_mass_A =   mp_estimated.*J_A./JSum;
        estimated_mass_productivity_A = estimated_mass_A.*mu_estimated_profile;
        estimated_specific_productivity_A = J_A./JSum.*mu_estimated_profile;
        
        r_estimated = mu_r_profile./mu_estimated_profile;
        rt_estimated = (1+JSum).*r_estimated/model_p.Phi_t;
        nu_t_si = f_substrate.*model_p.nu_max;
        BetaP_A_estimated = J_A.* r_estimated.*nu_t_si.*model_p.dm_A./model_p.lp_A./varying_omega_A(p); % Translation rate beta_p
        BetaP_R_estimated = JNr.* r_estimated.*nu_t_si.*model_p.dm_r./model_p.lp_r./model_p.Omega_r/model_p.N_r; % Translation rate beta_p
        BetaP_P_estimated = (JSum-JNr-J_A).* r_estimated.*nu_t_si.*model_p.dm_nr./model_p.lp_nr./model_p.Omega_nr/model_p.N_nr;
        
        % Each struct below gathers the results as f(s_i)-dependent lines for
        % each promoter and RBS strengths
        Struct_results_RBS.protein_chars{p,k}=proteinA_params;
        Struct_results_RBS.fraction_A{p,k}=estimated_fraction_A;
        Struct_results_RBS.fraction_R{p,k}=estimated_fraction_R;
        Struct_results_RBS.fraction_P{p,k}=estimated_fraction_P;
        Struct_results_RBS.mass_A{p,k}=estimated_mass_A;
        Struct_results_RBS.mass_productivity_A{p,k}=estimated_mass_productivity_A;
        Struct_results_RBS.specific_productivity_A{p,k}=estimated_specific_productivity_A;
        Struct_results_RBS.r{p,k} = r_estimated;
        Struct_results_RBS.rt{p,k} = rt_estimated;
        Struct_results_RBS.mu{p,k} = mu_estimated_profile;
        Struct_results_RBS.KC0_r{p,k} = KC0_r;
        Struct_results_RBS.KC0_nr{p,k} = KC0_nr;
        Struct_results_RBS.KC0_A{p,k} = KC0_A;
        Struct_results_RBS.BetaP_A{p,k} = BetaP_A_estimated;
        Struct_results_RBS.BetaP_P{p,k} = BetaP_P_estimated;
        Struct_results_RBS.BetaP_R{p,k} = BetaP_R_estimated;
        Struct_results_RBS.J_A{p,k} = J_A;
        Struct_results_RBS.flux_free_resources{p,k} = mu_r_profile;
    end
end

% The structs below gather the results as RBS-dependent lines for each
% substrate f(s_i) and each promoter strength
% We plot results for indices:
%

%%%%%%%%%%%%  PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures_plot_promoter=false;
figures_plot_RBS=true;
figures_plot_surfaces=true;

%%%%%% PLOTS FOR PROMOTER VARIATION:

%%%%%% PLOTS FOR RBS VARIATION:

f_substrate_round=round(f_substrate, 2);
% 
% f1 = figure(1);
% 
% substrate_choice_index=1; % correponds to f(s_i)=1
% RBS_strengths =check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index); %increasing RBS strengths. Log10 scale
% Promoter_strengths=varying_omega_A(1:4); %(1:4); %plotting_omega_A_indices; %increasing promoter strengths.
% [X_promoter1,Y_rbs1] = meshgrid(Promoter_strengths,RBS_strengths);
% Z_flux_free_resources1=zeros(length(Promoter_strengths),length(RBS_strengths));
% for i=1:length(Promoter_strengths)
%     for j=1:length(RBS_strengths)
%         Z_flux_free_resources1(i,j) = Struct_results_RBS.flux_free_resources{i,j}(substrate_choice_index);
%     end
% end
% 
% substrate_choice_index=6; % correponds to f(s_i)=1
% RBS_strengths =check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index); %increasing RBS strengths. Log10 scale
% Promoter_strengths=varying_omega_A; %(1:4); %plotting_omega_A_indices; %increasing promoter strengths.
% [X_promoter2,Y_rbs2] = meshgrid(Promoter_strengths,RBS_strengths);
% Z_flux_free_resources2=zeros(length(Promoter_strengths),length(RBS_strengths));
% for i=1:length(Promoter_strengths)
%     for j=1:length(RBS_strengths)
%         Z_flux_free_resources2(i,j) = Struct_results_RBS.flux_free_resources{i,j}(substrate_choice_index);
%     end
% end
% 
% substrate_choice_index=13; % correponds to f(s_i)=1
% RBS_strengths =check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index); %increasing RBS strengths. Log10 scale
% Promoter_strengths=varying_omega_A; %(1:4); %plotting_omega_A_indices; %increasing promoter strengths.
% [X_promoter3,Y_rbs3] = meshgrid(Promoter_strengths,RBS_strengths);
% Z_flux_free_resources3=zeros(length(Promoter_strengths),length(RBS_strengths));
% for i=1:length(Promoter_strengths)
%     for j=1:length(RBS_strengths)
%         Z_flux_free_resources3(i,j) = Struct_results_RBS.flux_free_resources{i,j}(substrate_choice_index);
%     end
% end
% 
% subplot(131)
% num_levels_contour=15; 
% [M_flux_free,c_flux_free]=contourf(...
%         X_promoter1,...
%         model_p.dm_A./Y_rbs1,...
%         Z_flux_free_resources1',...
%         num_levels_contour);
%         
% c_flux_free.LineWidth = 2;
% grid on, ylabel('$d_{m,A}/K_{C0}^A(s_i)  (\mathrm{molec}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex')
% xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
% hbc=colorbar;
% hbc.TickLabelInterpreter='latex';
% hbc.Title.String = '\mu r';
% hbc.Label.FontSize = 12;
% hbc.Title.FontSize = 12;
% substrate_choice_index=1;
% substrate_value=num2str(f_substrate_round(1,substrate_choice_index));
% s_title = strcat('f(s_i) =  ',substrate_value);
% title(s_title)
% 
% subplot(132)
% num_levels_contour=15;
% [M_flux_free,c_flux_free]=contourf(X_promoter2,model_p.dm_A./Y_rbs2,Z_flux_free_resources2',num_levels_contour);
% c_flux_free.LineWidth = 2;
% grid on, ylabel('$d_{m,A}/K_{C0}^A(s_i)  (\mathrm{molec}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex')
% xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
% hbc=colorbar;
% hbc.TickLabelInterpreter='latex';
% hbc.Title.String = '\mu r';
% hbc.Label.FontSize = 12;
% hbc.Title.FontSize = 12;
% substrate_choice_index=6;
% substrate_value=num2str(f_substrate_round(1,substrate_choice_index));
% s_title = strcat('f(s_i) =  ',substrate_value);
% title(s_title)
% 
% subplot(133)
% num_levels_contour=15;
% [M_flux_free,c_flux_free]=contourf(X_promoter3,model_p.dm_A./Y_rbs3,Z_flux_free_resources3',num_levels_contour);
% c_flux_free.LineWidth = 2;
% grid on, ylabel('$d_{m,A}/K_{C0}^A(s_i)  (\mathrm{molec}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex')
% xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
% hbc=colorbar;
% hbc.TickLabelInterpreter='latex';
% hbc.Title.String = '\mu r';
% hbc.Label.FontSize = 12;
% hbc.Title.FontSize = 12;
% substrate_choice_index=13;
% substrate_value=num2str(f_substrate_round(1,substrate_choice_index));
% s_title = strcat('f(s_i) =  ',substrate_value);
% title(s_title)

%%

ku_r = 130.74;
ku_nr = 3.11;
kb_r = 5.81;
kb_nr = 11.80;

KrC0 = kb_r / (ku_r + model_p.nu_max/model_p.l_e);
KnrC0 = kb_nr / (ku_nr + model_p.nu_max/model_p.l_e);

f2 = figure('Position',[500 500 500 400]);
clf(f2);

substrate_choice_index=1; % correponds to f(s_i)=1
RBS_strengths =check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index); %increasing RBS strengths. Log10 scale
Promoter_strengths=varying_omega_A; %(1:4); %plotting_omega_A_indices; %increasing promoter strengths.
[X_promoter1,Y_rbs1] = meshgrid(Promoter_strengths,RBS_strengths);
Z_flux_free_resources1=zeros(length(Promoter_strengths),length(RBS_strengths));
J_A = zeros(length(Promoter_strengths),length(RBS_strengths));
for i=1:length(Promoter_strengths)
    for j=1:length(RBS_strengths)
        Z_flux_free_resources1(i,j) = Struct_results_RBS.flux_free_resources{i,j}(substrate_choice_index);
        J_A(i,j) = Struct_results_RBS.J_A{i,j}(substrate_choice_index);
    end
end

kA = 1;
hA = model_p.dm_A ./ Y_rbs1;

% mur = 10*(0:0.1:1).*ones(121,1);
mur = [0.1 0.5 1 1 1 5 10 50 100 100 100].*ones(121,1);
% mur = [0.01 0.05 0.1 0.5 1 5 10 50 100 500 1000].*ones(121,1);

JA = kA./hA .* hA ./ (hA + mur);


num_levels_contour=15; 
[M_flux_free,c_flux_free]=contourf(...
        mur,...
        Y_rbs1,...
        JA,...
        num_levels_contour,'LineStyle','none');
        
c_flux_free.LineWidth = 0.5;
grid on, ylabel('$K_{C0}^A(s_i)  \quad [\mathrm{molec}^{-1}]$','Interpreter','latex')
xlabel('$log_{10}(\mu r) \quad [\mathrm{molec}\cdot \mathrm{min}^{-1}]$','Interpreter','latex');
hbc=colorbar;
title('$J_A(\mu,r) \quad [adim]$','Interpreter','latex')
ylim([0 0.3]);
set(gca, 'XScale', 'log')

caxis([min(min(JA)) max(max(JA))]);
hold on;
grid off;
    
% RBS ribosomal
plot(mur(1,:), KrC0*ones(size(X_promoter1(1,:))),'k--','LineWidth',1.5);
% RBS non-ribosomal
plot(mur(1,:), KnrC0*ones(size(X_promoter1(1,:))),'w--','LineWidth',1.5);

saveFig(f2,'fig_constantSubstrate_variableMuR',1,'/home/nobel/Sync/phd/paper/2020_Minimal_Bioarchive_Texte/images/');
