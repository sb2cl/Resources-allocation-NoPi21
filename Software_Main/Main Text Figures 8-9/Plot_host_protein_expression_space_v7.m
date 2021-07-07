%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state interactions host-protein considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp5.m for Nr=57
%   as analysed in >Analysis_resultys_estim_Nr57_v6.m and Plot_results_average_Nr57_v6.m 
%   for the host dynamics
%  Here we get the values of productivity as a function of promoter and RBS
%  strengths for the protein of interest (A) for varying substrate
%  availability
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
td_points=[30:5:100];
Bremer_exp_data.mu= log(2)./td_points;  %Row vector
b0=4.45687e-05;
b1=19.9843;
a1=34.9866;
gamma1= 0.679116;
gamma2=1.06627;
Bremer_exp_data.f_si= (b0 + b1*Bremer_exp_data.mu.^gamma1)./(1+a1*Bremer_exp_data.mu.^gamma2);

f_substrate = Bremer_exp_data.f_si; % ensure that goes from 0 to 1. Row vector

%%%%%%% General values of cell parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

model_p.l_e = 25; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.le_nr=model_p.lp_nr^0.097/0.0703;  %More precision 24.9874
model_p.le_r=model_p.lp_r^0.097/0.0703;  %More precision  23.7235
model_p.nu_max = 1260; % Maximum effective translation rate per ribosome (aa/min)
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57; %Number of protein types explaining 99% J's in Hausser data (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)		
model_p.Em_r = model_p.lp_r/ model_p.le_r*(1- ( model_p.lp_r/( model_p.lp_r+ model_p.le_r ) )^(model_p.lp_r/model_p.le_r)); 
model_p.Em_nr = model_p.lp_nr/ model_p.le_nr*(1- ( model_p.lp_nr/( model_p.lp_nr+ model_p.le_nr ) )^(model_p.lp_nr/model_p.le_nr)); 
model_p.WEm_r =  1 + 1/model_p.Em_r; %1.2891; % Average weight (1+1/Emr) for ribosomal protein-coding genes
model_p.WEm_nr =  1 + 1/model_p.Em_nr; % 1.1575; Average weight (1+1/Emp) for non-ribosomal protein-coding genes
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa)
Da = 1.6605e-24; % g/Da
model_p.ribosome_mass = 2.29*model_p.m_aa*model_p.N_r*model_p.lp_r; % mp_estimated includes the RNA component of ribosomes.

%%%%% Weighted average values of estimated parameters for the host

model_p.ku_r =  129.9034; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = 3.0921; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  5.5745; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  12.8570; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = 5.6509; % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = 0.02750; % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.Phi_m = 0.9043; % Average ratio between mature available and total number of ribosomes 

best_N57_params=[model_p.ku_r,model_p.ku_nr, model_p.kb_r, model_p.kb_nr, model_p.Omega_r, model_p.Omega_nr,model_p.Phi_m ];

model_p.Phi_s_b = 0.99*ones(1,length(f_substrate)); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_h_t = 0.83*ones(1,length(f_substrate)); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_r_t = 0.3*ones(1,length(f_substrate));  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_nr_t = 0.55*ones(1,length(f_substrate));  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_m = 0.91; % Initial estimation of the fraction of mature to inmature ribosomes

model_p.mu_estimated = Bremer_exp_data.mu; %Initial estimations for the function eval_host_protein_ss_v7.m . % Row vector

%%%%%%% Initial values for the promoter and RBS strengths of the exogenous
%%%%%%% protein A

proteinA='average'; %'galactosidase';

if strcmp(proteinA,'galactosidase') 
    %%%%%%% Initial values for the promoter and RBS strengths of the exogenous
    %%%%%%% protein A = beta-galactosidase. A large protein
    model_p.ku_A = 3.09; %Dissotiation rate RBS-ribosome  (1/min)            
    model_p.kb_A =  13.41; %Assotiation rate RBS-ribosome  (1/min/molec)
    model_p.Omega_A = 0.03; % Transcription rate   (1/min)           
    model_p.lp_A = 1023; % Assume protein A with mean characteristics
    model_p.dm_A = model_p.dm_nr;
    model_p.le_A= 27.86 ;   %lea=la^0.097/0.0703
    model_p.Em_A = model_p.lp_A/ model_p.le_A*(1- ( model_p.lp_A/( model_p.lp_A+ model_p.le_A ) )^(model_p.lp_A/model_p.le_A)); 
    model_p.WEm_A =  1 + 1/model_p.Em_A; %1.034; %Weight (1+1/EmA) EmA=23.0280
    proteinA_params =[model_p.ku_A, model_p.kb_A, model_p.Omega_A,model_p.lp_A,model_p.dm_A,model_p.le_A,model_p.WEm_A];
elseif strcmp(proteinA,'lactamase')   
    %%%%%%% Initial values for the promoter and RBS strengths of the exogenous
    %%%%%%% protein B = beta-lactamase. A small protein
    model_p.ku_A = 3.09; %Dissotiation rate RBS-ribosome  (1/min)            
    model_p.kb_A =  13.41; %Assotiation rate RBS-ribosome  (1/min/molec)
    model_p.Omega_A = 0.03; % Transcription rate   (1/min)           
    model_p.lp_A = 286; % 
    model_p.dm_A = model_p.dm_nr;
    model_p.le_A= 24.62 ;   %lea=la^0.097/0.0703
    model_p.Em_A = model_p.lp_A/ model_p.le_A*(1- ( model_p.lp_A/( model_p.lp_A+ model_p.le_A ) )^(model_p.lp_A/model_p.le_A)); 
    model_p.WEm_A =  1 + 1/model_p.Em_A; %1.1396; %Weight (1+1/EmA) EmA=7.1651
    proteinA_params =[model_p.ku_A, model_p.kb_A, model_p.Omega_A,model_p.lp_A,model_p.dm_A,model_p.le_A,model_p.WEm_A];
else
    % An average protein
    model_p.ku_A = 3.09; %Dissotiation rate RBS-ribosome  (1/min)            
    model_p.kb_A = 13.41; %Assotiation rate RBS-ribosome  (1/min/molec)
    model_p.Omega_A = 0.03; % Transcription rate   (1/min)           
    model_p.lp_A = model_p.lp_nr; % Assume protein A with mean characteristics
    model_p.le_A= model_p.le_nr;   %lea=la^0.097/0.0703
    model_p.dm_A = model_p.dm_nr;
    model_p.Em_A =  model_p.Em_nr;
    model_p.WEm_A =  model_p.WEm_nr; %Weight (1+1/EmA) for average non-ribosomal protein-coding genes. 
    proteinA_params =[model_p.ku_A, model_p.kb_A, model_p.Omega_A,model_p.lp_A,model_p.dm_A,model_p.le_A,model_p.WEm_A];    
end

%%% Cell protein mass and dry weight

model_p.model_mass=1;

%%%%% EVALUATION FOR AVERAGE RBS AND VARYING PROMOTER:

varying_omega_A=[0,1,10,25, 50, 100, 250, 500,750, 1000,1500, 2000,2500, 5000,7500,10000, 12500, 15000,20000, 25000,35000,50000,65000,80000, 120000]*model_p.Omega_A;

for k=1:length(varying_omega_A)
    
    proteinA_params =[model_p.ku_A, model_p.kb_A, varying_omega_A(k),model_p.lp_A,model_p.dm_A,model_p.le_A,model_p.WEm_A];

    [ke_A,KC0_r,KC0_nr, KC0_A,mu_estimated_profile,estimated_mass_h,mu_r_profile, JNr, JNnr, JNA, JWSum]   = eval_host_protein_ss_v7(proteinA_params, f_substrate); 
    
     JSum = JNr + JNnr + JNA;
     
    % Strain mass:
     estimated_mass_s = JSum./(JNr + JNnr).*estimated_mass_h;
     
     % Estimated fractions relative to the strain  mass:
    estimated_fraction_A_s =   JNA./JSum;
    estimated_fraction_R_s = JNr./JSum;
    estimated_fraction_P_s = 1-estimated_fraction_A_s-estimated_fraction_R_s;
    
    % Estimated masses:
    estimated_mass_A =   estimated_mass_s.*estimated_fraction_A_s;
    estimated_mass_R =   estimated_mass_s.*estimated_fraction_R_s;
    estimated_mass_P =   estimated_mass_s.*estimated_fraction_P_s;
    
    % Estimated mass of fraction Q:
    if model_p.model_mass==1 % this is the mass of the native host
        mc0=123.022e-15; %(g)  
        gamma=68.5547; % (min)
        estimated_mass_cDW = mc0*exp(gamma*mu_estimated_profile) + estimated_mass_A; 
    else
       % We do not use anothe one for the cell DW, as they tend to give
       % smaller values than the protein mass for very low growth rate with
       % the experimental data we had to fit the models. See
       % Main_fit_Mass_h_mu_2p.m and Main_fit_ProteinMass_h_mu.m
        mc0=123.022e-15; %(g)  
        gamma=68.5547; % (min)
        estimated_mass_cDW = mc0*exp(gamma*mu_estimated_profile) + estimated_mass_A; 
    end
    
    % Estimated fractions relative to the cDW mass:
    estimated_fraction_A_cDW =  estimated_mass_A./estimated_mass_cDW;
    estimated_fraction_R_cDW = estimated_mass_R./estimated_mass_cDW;
    estimated_fraction_P_cDW = estimated_mass_P./estimated_mass_cDW;
    estimated_fraction_Q_cDW = 1-estimated_fraction_A_cDW-estimated_fraction_R_cDW-estimated_fraction_P_cDW;
    
    % Estimated synthesis rates:
     estimated_mass_productivity_A = estimated_mass_A.*mu_estimated_profile;
     estimated_specific_productivity_A = estimated_mass_productivity_A./estimated_mass_cDW;  %per cell dry weight
    
    % Estimated number of free and total ribosomes 
    r_estimated = mu_r_profile./mu_estimated_profile; %r_estimated = estimated_flux_mur/estimated_mu
    rT_estimated = (1+JWSum).*r_estimated/model_p.Phi_m; %Estimated total number of ribosomes: r_a = Phi_m*r_T -> r_T=r_a/Phi_m = r*(1+JWSum)/Phi_m
    
    % Store data:
    Struct_results_promoter.protein_chars{k}=proteinA_params;
    Struct_results_promoter.fraction_A{k}=estimated_fraction_A_cDW;
    Struct_results_promoter.fraction_R{k}=estimated_fraction_R_cDW;
    Struct_results_promoter.fraction_P{k}=estimated_fraction_P_cDW;
    Struct_results_promoter.fraction_Q{k}=estimated_fraction_Q_cDW;
    Struct_results_promoter.mass_A{k}=estimated_mass_A;
    Struct_results_promoter.mass_productivity_A{k}=estimated_mass_productivity_A;
    Struct_results_promoter.specific_productivity_A{k}=estimated_specific_productivity_A;
    Struct_results_promoter.r{k} = r_estimated;
    Struct_results_promoter.rt{k} = rT_estimated;
    Struct_results_promoter.mu{k} = mu_estimated_profile;
end

%%%%% EVALUATION FOR AVERAGE  PROMOTER  AND VARYING RBS:
% Recall average values model_p.ku_A = 3.09; %Dissotiation rate RBS-ribosome  (1/min)            
% model_p.kb_A =  13.41; %Assotiation rate RBS-ribosome  (1/min/molec)
range_Ku=[3,135]; 
range_Kb=[3,15];
n_levels=10;
discrete_step_Ku=round((max(range_Ku)-min(range_Ku))/n_levels,0);
discrete_step_Kb=round((max(range_Kb)-min(range_Kb))/n_levels,1);
dFF = fullfact([n_levels+1, n_levels+1]); % Full factorial experiment
varying_matrix_Ku_Kb_A = [3+(dFF(:,1)-1)*discrete_step_Ku, 3+(dFF(:,2)-1)*discrete_step_Kb];
varying_ratio_Kb_over_Ku_A= varying_matrix_Ku_Kb_A(:,2)./(varying_matrix_Ku_Kb_A(:,1)+model_p.nu_max/model_p.le_A); %RBS strength for substrate saturation
[sorted_Kb_over_Ku_A,index_sorted_Kb_over_Ku_A]=sort(varying_ratio_Kb_over_Ku_A,'ascend'); %sorted RBS strengths
sorted_varying_matrix_Ku_Kb_A = varying_matrix_Ku_Kb_A(index_sorted_Kb_over_Ku_A,:); % sorted matrix of (Ku,Kb) pairs according to ascending RBS strength for saturated substrate
%%%%%%% Not needed to check, for there is a monotonous relationship, but gives us the
%%%%%%% set of all RBS strengths for all values of the substrate concentration function
%%%%%%% f(s_i) needed for plotting 
 check_varying_ratio_Kb_over_Ku_A=[];
 for k=1:length(f_substrate)
     check_varying_ratio_Kb_over_Ku_A(:,k)= sorted_varying_matrix_Ku_Kb_A(:,2)./(sorted_varying_matrix_Ku_Kb_A(:,1)+model_p.nu_max/model_p.le_A*f_substrate(k)); 
 end
% check_varying_ratio_Kb_over_Ku_A
%%%%%%%%%%%%%%%%%%%%%%%%%

%varying_omega_A=[0,1,10,25, 50, 250, 500,750, 1000,1500, 2000,2500, 5000,7500,10000, 12500, 15000,20000, 25000,35000,50000,65000,80000, 120000]*model_p.Omega_A;


for p=1:length(varying_omega_A)
    for k=1:size(sorted_varying_matrix_Ku_Kb_A,1)
    % Recall average model_p.Omega_A = 0.03; % Transcription rate   (1/min)     
    
    proteinA_params =[sorted_varying_matrix_Ku_Kb_A(k,1),sorted_varying_matrix_Ku_Kb_A(k,2),varying_omega_A(p),model_p.lp_A,model_p.dm_A,model_p.le_A,model_p.WEm_A];

    [ke_A,KC0_r,KC0_nr, KC0_A,mu_estimated_profile,estimated_mass_h,mu_r_profile, JNr, JNnr, JNA, JWSum]   = eval_host_protein_ss_v7(proteinA_params, f_substrate); 
    
   JSum = JNr + JNnr + JNA;
     
    % Strain mass:
     estimated_mass_s = JSum./(JNr + JNnr).*estimated_mass_h;
     
     % Estimated fractions relative to the strain  mass:
    estimated_fraction_A_s =   JNA./JSum;
    estimated_fraction_R_s = JNr./JSum;
    estimated_fraction_P_s = 1-estimated_fraction_A_s-estimated_fraction_R_s;
    
    % Estimated masses:
    estimated_mass_A =   estimated_mass_s.*estimated_fraction_A_s;
    estimated_mass_R =   estimated_mass_s.*estimated_fraction_R_s;
    estimated_mass_P =   estimated_mass_s.*estimated_fraction_P_s;
    
    % Estimated mass of fraction Q:
    if model_p.model_mass==1 % this is the mass of the native host
        mc0=123.022e-15; %(g)  
        gamma=68.5547; % (min)
        estimated_mass_cDW = mc0*exp(gamma*mu_estimated_profile) + estimated_mass_A; 
    else
       % We do not use anothe one for the cell DW, as they tend to give
       % smaller values than the protein mass for very low growth rate with
       % the experimental data we had to fit the models. See
       % Main_fit_Mass_h_mu_2p.m and Main_fit_ProteinMass_h_mu.m
        mc0=123.022e-15; %(g)  
        gamma=68.5547; % (min)
        estimated_mass_cDW = mc0*exp(gamma*mu_estimated_profile) + estimated_mass_A; 
    end
    
    % Estimated fractions relative to the cDW mass:
    estimated_fraction_A_cDW =  estimated_mass_A./estimated_mass_cDW;
    estimated_fraction_R_cDW = estimated_mass_R./estimated_mass_cDW;
    estimated_fraction_P_cDW = estimated_mass_P./estimated_mass_cDW;
    estimated_fraction_Q_cDW = 1-estimated_fraction_A_cDW-estimated_fraction_R_cDW-estimated_fraction_P_cDW;
    
    % Estimated synthesis rates:
     estimated_mass_productivity_A = estimated_mass_A.*mu_estimated_profile;
     estimated_specific_productivity_A = estimated_mass_productivity_A./estimated_mass_cDW;  %per cell dry weight
    
    % Estimated number of free and total ribosomes 
    r_estimated = mu_r_profile./mu_estimated_profile; %r_estimated = estimated_flux_mur/estimated_mu
    rT_estimated = (1+JWSum).*r_estimated/model_p.Phi_m; %Estimated total number of ribosomes: r_a = Phi_m*r_T -> r_T=r_a/Phi_m = r*(1+JWSum)/Phi_m
    
    % Each struct below gathers the results as f(s_i)-dependent lines for
    % each promoter and RBS strengths
    Struct_results_RBS.protein_chars{p,k}=proteinA_params;
    Struct_results_RBS.fraction_A{p,k}=estimated_fraction_A_cDW;
    Struct_results_RBS.fraction_R{p,k}=estimated_fraction_R_cDW;
    Struct_results_RBS.fraction_P{p,k}=estimated_fraction_P_cDW;
    Struct_results_RBS.fraction_Q{p,k}=estimated_fraction_Q_cDW;

    Struct_results_RBS.mass_A{p,k}=estimated_mass_A;
    Struct_results_RBS.mass_productivity_A{p,k}=estimated_mass_productivity_A;
    Struct_results_RBS.specific_productivity_A{p,k}=estimated_specific_productivity_A;
    Struct_results_RBS.r{p,k} = r_estimated;
    Struct_results_RBS.rt{p,k} = rT_estimated;
    Struct_results_RBS.mu{p,k} = mu_estimated_profile;
    
    Struct_results_RBS.KC0_r{p,k} = KC0_r;
    Struct_results_RBS.KC0_nr{p,k} = KC0_nr;
    Struct_results_RBS.KC0_A{p,k} = KC0_A;
   
    end
end

% The structs below gather the results as RBS-dependent lines for each
% substrate f(s_i) and each promoter strength
% Recall 
% varying_omega_A=[0,1,10,25, 50, 250, 500,750, 1000,1500, 2000,2500, 5000,7500,10000, 12500, 15000,20000, 25000,35000,50000,65000,80000, 120000]*model_p.Omega_A;

% We plot results for indices:
%
plotting_omega_A_indices = [10,20]; 

for p=1:length(plotting_omega_A_indices)
    for m=1:length(f_substrate)
         lines_fraction_A_si.dades{p,m}=[];
         lines_fraction_P_si.dades{p,m}=[];
         lines_fraction_R_si.dades{p,m}=[];
         lines_fraction_Q_si.dades{p,m}=[];
         lines_mu_si.dades{p,m}=[];
         lines_mass_productivity_A_si.dades{p,m}=[];
         lines_specific_productivity_A_si.dades{p,m}=[];
         for n=1:size(check_varying_ratio_Kb_over_Ku_A,1)
           lines_fraction_A_si.dades{p,m}=[lines_fraction_A_si.dades{p,m},Struct_results_RBS.fraction_A{plotting_omega_A_indices(p),n}(m)];
           lines_fraction_P_si.dades{p,m}=[lines_fraction_P_si.dades{p,m},Struct_results_RBS.fraction_P{plotting_omega_A_indices(p),n}(m)];
           lines_fraction_R_si.dades{p,m}=[lines_fraction_R_si.dades{p,m},Struct_results_RBS.fraction_R{plotting_omega_A_indices(p),n}(m)];
           lines_fraction_Q_si.dades{p,m}=[lines_fraction_Q_si.dades{p,m},Struct_results_RBS.fraction_Q{plotting_omega_A_indices(p),n}(m)];
           lines_mu_si.dades{p,m}=[lines_mu_si.dades{p,m},Struct_results_RBS.mu{plotting_omega_A_indices(p),n}(m)];
           lines_mass_productivity_A_si.dades{p,m}=[ lines_mass_productivity_A_si.dades{p,m},Struct_results_RBS.mass_productivity_A{plotting_omega_A_indices(p),n}(m)];
           lines_specific_productivity_A_si.dades{p,m}=[ lines_specific_productivity_A_si.dades{p,m},Struct_results_RBS.specific_productivity_A{plotting_omega_A_indices(p),n}(m)];
         end
    end
end

%%%%%%%%%%%%  PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

figures_plot_promoter=true;
figures_plot_RBS=true;
figures_plot_surfaces=true;

%%%%%% PLOTS FOR PROMOTER VARIATION:
if figures_plot_promoter==true
ColorMap = jet(length(f_substrate));
Linewidth = 3;
Fontsize = 16;
   
     f1= figure(1);
         subplot(131)
         for m=1:length(f_substrate)
             line_fraction_A_si=[];
             line_fraction_P_si=[];
             line_fraction_R_si=[];
             line_fraction_Q_si=[];
             for n=1:length(varying_omega_A)
               line_fraction_A_si=[line_fraction_A_si,Struct_results_promoter.fraction_A{n}(m)];
               line_fraction_P_si=[line_fraction_P_si,Struct_results_promoter.fraction_P{n}(m)];
               line_fraction_R_si=[line_fraction_R_si,Struct_results_promoter.fraction_R{n}(m)];
               line_fraction_Q_si=[line_fraction_Q_si,Struct_results_promoter.fraction_Q{n}(m)];
             end
             if m==length(f_substrate)
             plot(log10(varying_omega_A),  line_fraction_A_si,'Color',ColorMap(m,:),'MarkerSize',8,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(varying_omega_A),  line_fraction_R_si+line_fraction_P_si+line_fraction_Q_si,'Color',ColorMap(m,:),'MarkerSize',8,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','Diamond', 'LineStyle','-','LineWidth',2);
             else
                 plot(log10(varying_omega_A),  line_fraction_A_si,'Color',ColorMap(m,:),'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(varying_omega_A),  line_fraction_R_si+line_fraction_P_si+line_fraction_Q_si,'Color',ColorMap(m,:),'LineStyle','-','LineWidth',2);
             end      
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Phi_A\, (\mathrm{O}), \Phi_{NR} + \Phi_{R} + \Phi_{Q} (\diamondsuit)$','FontSize',Fontsize,'Interpreter','latex');
         hold off
          ax = gca;
         ax.FontSize = Fontsize;
         subplot(132)
         for m=1:length(f_substrate)
             line_mu_si=[];
             for n=1:length(varying_omega_A)
               line_mu_si=[line_mu_si,Struct_results_promoter.mu{n}(m)];
             end
             plot(log10(varying_omega_A), line_mu_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\mu\,  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
          ax = gca;
         ax.FontSize = Fontsize;
          subplot(133)
         for m=1:length(f_substrate)
            % line_mass_productivity_A_si=[];
             line_specific_productivity_A_si=[];
             for n=1:length(varying_omega_A)
               % line_mass_productivity_A_si=[ line_mass_productivity_A_si,Struct_results_promoter.mass_productivity_A{n}(m)];
                line_specific_productivity_A_si=[ line_specific_productivity_A_si,Struct_results_promoter.specific_productivity_A{n}(m)];
             end
             %plot(log10(varying_omega_A),  line_mass_productivity_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
             plot(log10(varying_omega_A), line_specific_productivity_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\pi_A\,  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{cDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
          ax = gca;
         ax.FontSize = Fontsize;
         colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
           hbc.Title.Interpreter = 'latex';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = '$f(s_i)$';
         hbc.Label.FontSize = 14;
         hbc.Title.FontSize = 12;
          
         %exportgraphics(f1,'./promoter_swap.png','Resolution',300)
          
         
           f2= figure(2);
         for m=1:length(f_substrate)
             line_fraction_A_si=[];
             line_fraction_P_si=[];
             line_fraction_R_si=[];
              line_fraction_Q_si=[];
             for n=1:length(varying_omega_A)
               line_fraction_A_si=[line_fraction_A_si,Struct_results_promoter.fraction_A{n}(m)];
               line_fraction_P_si=[line_fraction_P_si,Struct_results_promoter.fraction_P{n}(m)];
               line_fraction_R_si=[line_fraction_R_si,Struct_results_promoter.fraction_R{n}(m)];
               line_fraction_Q_si=[line_fraction_Q_si,Struct_results_promoter.fraction_Q{n}(m)];
             end
             if m==length(f_substrate)
             plot(log10(varying_omega_A),  line_fraction_A_si,'Color',ColorMap(m,:),'MarkerSize',12,...
                  'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(varying_omega_A),  line_fraction_R_si+line_fraction_P_si+line_fraction_Q_si,'Color',ColorMap(m,:),'MarkerSize',12,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','Diamond', 'LineStyle','-','LineWidth',2);
             else
                 plot(log10(varying_omega_A),  line_fraction_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(varying_omega_A),  line_fraction_R_si+line_fraction_P_si+line_fraction_Q_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
             end
         end
         ax = gca;
         ax.FontSize = Fontsize;
         grid on, xlabel('$\log_{10}{N_A\omega_A}\, (\mathrm{mRNA}\cdot\mathrm{min}^{-1})$','FontSize',Fontsize+8,'Interpreter','latex')
         ylabel('$\Phi_A\, (\mathrm{O}), \Phi_{NR} + \Phi_{R}+ \Phi_{Q} (\diamondsuit)$','FontSize',Fontsize+8,'Interpreter','latex');
         hold off
         %exportgraphics(f2,'../images_def/promoter_swap_subplot_fractions.png','Resolution',300)

         
          f3= figure(3);
         subplot(121)
         for m=1:length(f_substrate)
             line_mu_si=[];
             for n=1:length(varying_omega_A)
               line_mu_si=[line_mu_si,Struct_results_promoter.mu{n}(m)];
             end
             plot(log10(varying_omega_A), line_mu_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}\, (\mathrm{mRNA}\cdot\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\mu\,  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         axis([-2 4 0 0.03])
          ax = gca;
         ax.FontSize = Fontsize;
          colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
           hbc.Title.Interpreter = 'latex';
          hbc.Location='west';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = '$f(s_i)$';
         hbc.Label.FontSize = 16;
         hbc.Title.FontSize = 16;
          subplot(122)
         for m=1:length(f_substrate)
             line_specific_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_specific_productivity_A_si=[ line_specific_productivity_A_si,Struct_results_promoter.specific_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  line_specific_productivity_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}\, (\mathrm{mRNA}\cdot\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\pi_A  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{gCDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
          ax = gca;
         ax.FontSize = Fontsize;
         ax.YAxisLocation='right';
%          colormap(ColorMap);
%          f_substrate_round=round(f_substrate, 2);
%          hbc=colorbar('Direction','reverse');
%          hbc.TickLabelsMode='manual';
%          hbc.TickLabelInterpreter='latex';
%          hbc.TickLabels=f_substrate_round;
%          hbc.Title.String = 'f(s_i)';
%          hbc.Label.FontSize = 16;
%          hbc.Title.FontSize = 14;
          
          %exportgraphics(f3,'../images_def/promoter_swap_main_plot.png','Resolution',300)
         
         
         
          f4= figure(4);
           subplot(121)
         for m=1:length(f_substrate)
             line_mass_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_mass_productivity_A_si=[ line_mass_productivity_A_si,Struct_results_promoter.mass_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  1e15*line_mass_productivity_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{\omega_A}\, (\mathrm{mRNA}\cdot\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Pi_A  (fg\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
          ax = gca;
         ax.FontSize = Fontsize;
          subplot(122)
         for m=1:length(f_substrate)
             line_specific_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_specific_productivity_A_si=[ line_specific_productivity_A_si,Struct_results_promoter.specific_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  line_specific_productivity_A_si,'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{N_A\omega_A}\, (\mathrm{mRNA}\cdot\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\pi_A  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{gCDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
          ax = gca;
         ax.FontSize = Fontsize;
         colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
           hbc.Title.Interpreter = 'latex';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = '$f(s_i)$';
         hbc.Label.FontSize = 14;
         hbc.Title.FontSize = 12;
         
         %exportgraphics(f4,'../images_def/promoter_swap_productivities.png','Resolution',300)
end

%%%%%% PLOTS FOR RBS VARIATION:
if figures_plot_RBS==true
ColorMap = jet(length(f_substrate));
Linewidth = 3;
Fontsize = 16;
          
% f5= figure(5);
% row = 1;
% for p=1:length(plotting_omega_A_indices)*3
%     subplot(length(plotting_omega_A_indices),3,p)
%      if mod(p,3)==1
%          for m=1:length(f_substrate)
%              plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_fraction_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
%               hold on
%               plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  lines_fraction_P_si.dades{row,m}+lines_fraction_R_si.dades{row,m}+lines_fraction_Q_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
%          end
%           grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
%          ylabel('$\Phi_A (\mathrm{O}), \Phi_{NR} + \Phi_{R} + \Phi_{Q} (\diamondsuit)$','FontSize',Fontsize,'Interpreter','latex');
%          hold off
%          ax = gca;
%          ax.FontSize = Fontsize;
%      elseif mod(p,3)==2
%           for m=1:length(f_substrate)
%                 plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)), lines_mu_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
%                  hold on
%          end
%          grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
%          %axis([-1.8,-0.2,0,0.03])
%          ylabel('$\mu  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
%          omega_a_value=num2str(varying_omega_A(1,plotting_omega_A_indices(row)));
%          s_title = strcat('N_A \omega_A = ',omega_a_value);
%          title(s_title)
%          hold off
%          ax = gca;
%          ax.FontSize = Fontsize;
%      else
%         for m=1:length(f_substrate)
%             plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_specific_productivity_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
%             %plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  1e15*lines_mass_productivity_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
%              hold on
%          end
%          grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
%          %ylabel('$\Pi_A  (fg\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
%          ylabel('$\pi_A  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{gCDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
%          % axis([-2,-0.25,0,0.6])
%          hold off
%          ax = gca;
%          ax.FontSize = Fontsize;
%           colormap(ColorMap);
%          f_substrate_round=round(f_substrate, 2);
%          hbc=colorbar('Direction','reverse');
%          hbc.TickLabelsMode='manual';
%            hbc.Title.Interpreter = 'latex';
%          hbc.TickLabelInterpreter='latex';
%          hbc.TickLabels=f_substrate_round;
%          hbc.Title.String = '$f(s_i)$';
%          hbc.Label.FontSize = 2;
%          hbc.Title.FontSize = 16;
%           row = row+1;
%      end    
% end
% %exportgraphics(f5,'../images_def/RBS_swap_original.png','Resolution',300)

f5= figure(5);
row = 1;
for p=1:length(plotting_omega_A_indices)*2
    subplot(length(plotting_omega_A_indices),2,p)
     if mod(p,2)==1
         for m=1:length(f_substrate)
                plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)), lines_mu_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
                 hold on
         end
         grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         %axis([-1.8,-0.2,0,0.03])
         ylabel('$\mu  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         omega_a_value=num2str(varying_omega_A(1,plotting_omega_A_indices(row)));
         s_title = strcat('$N_A \omega_A =$',omega_a_value,'$\,\mathrm{mRNA}\cdot\mathrm{min}^{-1}$');
         title(s_title,'FontSize',Fontsize,'Interpreter','latex')
         hold off
         ax = gca;
         ax.FontSize = Fontsize;
     else
        for m=1:length(f_substrate)
            plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_specific_productivity_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
            %plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  1e15*lines_mass_productivity_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
             hold on
         end
         grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         %ylabel('$\Pi_A  (fg\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         ylabel('$\pi_A  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{gCDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         % axis([-2,-0.25,0,0.6])
         hold off
         ax = gca;
         ax.FontSize = Fontsize;
          colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
           hbc.Title.Interpreter = 'latex';
         hbc.TickLabelInterpreter='latex';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = '$f(s_i)$';
         hbc.Label.FontSize = 2;
         hbc.Title.FontSize = 16;
          row = row+1;
     end    
end
%exportgraphics(f5,'../images_def/RBS_swap_original.png','Resolution',300)

f6= figure(6);
row = 1;
for p=1:length(plotting_omega_A_indices)*2
    subplot(length(plotting_omega_A_indices),2,p)
       if mod(p,2)==1
          for m=1:length(f_substrate)
                plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)), lines_mu_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
                 hold on
         end
         grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)}\, (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         %axis([-1.8,-0.2,0,0.03])
         ylabel('$\mu\,  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         omega_a_value=num2str(varying_omega_A(1,plotting_omega_A_indices(row)));
         s_title = strcat('N_A \omega_A = ',omega_a_value);
        ax = gca;
        t= title(['$N_A \omega_A =$ ',omega_a_value,'$\,\mathrm{mRNA}\cdot\mathrm{min}^{-1}$'],'FontSize',Fontsize,'Interpreter','latex');
        t.FontWeight='normal';
         hold off
         ax.FontSize = Fontsize;
         colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
          hbc.Location='west';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = '$f(s_i)$';
         hbc.Label.FontSize = 16;
         hbc.Title.FontSize = 14;
     else
        for m=1:length(f_substrate)
            plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  1e15*lines_mass_productivity_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
             hold on
         end
         grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)}\, (\mathrm{molec}^{-1})$','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Pi_A\,  (fg\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
        %  axis([-2,-0.25,0,6e-16])
         hold off
         ax = gca;
         ax.FontSize = Fontsize;
          row = row+1;
     end    
end
%exportgraphics(f6,'../images_def/RBS_swap_main_plot.png','Resolution',300)

f7= figure(7);
row = 1;
        for m=1:length(f_substrate)
               if m==length(f_substrate)
             plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_fraction_A_si.dades{row,m},'Color','k','MarkerSize',12,...
                  'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  lines_fraction_P_si.dades{row,m}+lines_fraction_R_si.dades{row,m}+lines_fraction_Q_si.dades{row,m},'Color',ColorMap(m,:),'MarkerSize',12,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','Diamond', 'LineStyle','-','LineWidth',2);
               else
                   plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_fraction_A_si.dades{row,m},'Color',ColorMap(m,:)', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  lines_fraction_P_si.dades{row,m}+lines_fraction_R_si.dades{row,m}+lines_fraction_Q_si.dades{row,m},'Color',ColorMap(m,:),'LineStyle','-','LineWidth',2);
               end
               
         end
          grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize+8,'Interpreter','latex')
         ylabel('$\Phi_A (\mathrm{O}), \Phi_{NR} + \Phi_{R}  + \Phi_{Q} (\diamondsuit)$','FontSize',Fontsize+8,'Interpreter','latex');
         hold off
         ax = gca;
         ax.FontSize = Fontsize+8;
    %exportgraphics(f7,'../images_def/RBS_swap_subplot_fraction1.png','Resolution',300)
    
    

f8= figure(8);
row = 2;
        for m=1:length(f_substrate)
             if m==length(f_substrate)
             plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_fraction_A_si.dades{row,m},'Color','k','MarkerSize',12,...
                  'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  lines_fraction_P_si.dades{row,m}+lines_fraction_R_si.dades{row,m}+lines_fraction_Q_si.dades{row,m},'Color',ColorMap(m,:),'MarkerSize',12,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','Diamond', 'LineStyle','-','LineWidth',2);
             else
                plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),   lines_fraction_A_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(check_varying_ratio_Kb_over_Ku_A(:,m)),  lines_fraction_P_si.dades{row,m}+lines_fraction_R_si.dades{row,m}+lines_fraction_Q_si.dades{row,m},'Color',ColorMap(m,:), 'LineStyle','-','LineWidth',2); 
             end
         end
          grid on, xlabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',Fontsize+8,'Interpreter','latex')
         ylabel('$\Phi_A (\mathrm{O}), \Phi_{NR} + \Phi_{R} + \Phi_{Q} (\diamondsuit)$','FontSize',Fontsize+8,'Interpreter','latex');
         hold off
         ax = gca;
         ax.FontSize = Fontsize+8;
    %exportgraphics(f8,'../images_def/RBS_swap_subplot_fraction2.png','Resolution',300)

    

end

%%%%% PLOTS level lines:

if figures_plot_surfaces==true
ColorMap = jet(100);
Linewidth = 3;
Fontsize = 8;
          
substrate_choice_index=1; % correponds to f(s_i)=1
% RBS_strengths =log10(check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index)); %increasing RBS strengths. Log10 scale
% Promoter_strengths=log10(varying_omega_A); %plotting_omega_A_indices;
% %increasing promoter strengths.  Log10 scale

RBS_strengths =check_varying_ratio_Kb_over_Ku_A(:,substrate_choice_index); %increasing RBS strengths. Log10 scale
Promoter_strengths=varying_omega_A; %(1:4); %plotting_omega_A_indices; %increasing promoter strengths.

[X_promoter,Y_rbs] = meshgrid(Promoter_strengths,RBS_strengths);
% Now for each promoter-rbs strength pair we get the corresponding growth
% rate and productivities
Z_mu=zeros(length(Promoter_strengths),length(RBS_strengths));
for i=1:length(Promoter_strengths)
    for j=1:length(RBS_strengths)
            Z_mu(i,j) = Struct_results_RBS.mu{i,j}(substrate_choice_index);
    end
end

Z_mass_productivity_A=zeros(length(Promoter_strengths),length(RBS_strengths));
for i=1:length(Promoter_strengths)
    for j=1:length(RBS_strengths)
            Z_mass_productivity_A(i,j) = Struct_results_RBS.mass_productivity_A{i,j}(substrate_choice_index);
    end
end

Z_specific_productivity_A=zeros(length(Promoter_strengths),length(RBS_strengths));
for i=1:length(Promoter_strengths)
    for j=1:length(RBS_strengths)
            Z_specific_productivity_A(i,j) = Struct_results_RBS.specific_productivity_A{i,j}(substrate_choice_index);
    end
end


f9= figure(9);
subplot(121)
num_levels_contour=15;
% Z = sin(X)+cos(Y);
[M_mu,c_mu]=contour(X_promoter,Y_rbs,Z_mu',num_levels_contour);
c_mu.LineWidth = 2;
grid on, ylabel('$\log_{10}{K_{C0}^A(s_i)}  (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = '$\mu$';
 hbc.Label.FontSize = 8;
 hbc.Title.FontSize = 12;

subplot(122)
num_levels_contour=15;
[M_td,c_td]=contour(X_promoter,Y_rbs,log(2)./Z_mu',num_levels_contour);
c_td.LineWidth = 2;
grid on, ylabel('\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = '$t_d\,(\mathrm{min}))';
 hbc.Label.FontSize = 8;
 hbc.Title.FontSize = 12;
 
 f10= figure(10);
subplot(121)
num_levels_contour=15;
[M_mass_productivity_A,c_mass_productivity_A]=contourf(X_promoter,Y_rbs,Z_mass_productivity_A',num_levels_contour);
c_mass_productivity_A.LineWidth = 2;
grid on, ylabel('$\log_{10}{K_{C0}^A(s_i)}  (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = '$\Pi_A$';
 hbc.Label.FontSize = 12;
 hbc.Title.FontSize = 12;

subplot(122)
num_levels_contour=15;
[M_specific_productivity_A,c_specific_productivity_A]=contourf(X_promoter,Y_rbs,Z_specific_productivity_A',num_levels_contour);
hold on
plot(model_p.Omega_r*model_p.N_r,Struct_results_RBS.KC0_r{1,1}(substrate_choice_index) ,'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','d')                 
plot(model_p.Omega_nr*model_p.N_nr,Struct_results_RBS.KC0_nr{1,1}(substrate_choice_index),'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','s')  
c_specific_productivity_A.LineWidth = 2;
grid on, ylabel('$\log_{10}{K_{C0}^A(s_i)} (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = '$\pi_A\,  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{gCDW}^{-1})$';
 hbc.Label.FontSize = 12;
 hbc.Title.FontSize = 12;


 f11= figure(11);
 subplot(121)
num_levels_contour=15;
[M_mu,c_mu]=contourf(X_promoter,Y_rbs,log(2)./Z_mu',num_levels_contour);
hold on
plot(model_p.Omega_r*model_p.N_r,Struct_results_RBS.KC0_r{1,1}(substrate_choice_index) ,'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','d')                 
plot(model_p.Omega_nr*model_p.N_nr,Struct_results_RBS.KC0_nr{1,1}(substrate_choice_index),'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','s')  
c_mu.LineWidth = 2;
grid on, ylabel('$K_{C0}^A(s_i) (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = 't_d';
 hbc.Label.FontSize = 8;
 hbc.Title.FontSize = 12;
 hold off
 
subplot(122)
num_levels_contour=15;
[M_mass_productivity_A,c_mass_productivity_A]=contourf(X_promoter,Y_rbs,Z_mass_productivity_A',num_levels_contour);
hold on
plot(model_p.Omega_r*model_p.N_r,Struct_results_RBS.KC0_r{1,1}(substrate_choice_index) ,'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','d')                 
plot(model_p.Omega_nr*model_p.N_nr,Struct_results_RBS.KC0_nr{1,1}(substrate_choice_index),'Color','k','MarkerSize',18,'MarkerFaceColor','k', 'MarkerEdgeColor',[0 0 0], 'Marker','s')  
c_mass_productivity_A.LineWidth = 2;
grid on, ylabel('$K_{C0}^A(s_i)  (\mathrm{molec}^{-1})$','FontSize',14,'Interpreter','latex')
xlabel('$N_a\omega_A  (\mathrm{mRNA}\cdot \mathrm{min}^{-1})$','FontSize',14,'Interpreter','latex');
hbc=colorbar;
 hbc.TickLabelInterpreter='latex';
  hbc.Title.Interpreter = 'latex';
 hbc.Title.String = '$\Pi_A\,  (fg\cdot \mathrm{min}^{-1})$';
 hbc.Label.FontSize = 8;
 hbc.Title.FontSize = 12;
hold off

%exportgraphics(f6,'./expression_space_swap_fsimax.png','Resolution',300)



end




