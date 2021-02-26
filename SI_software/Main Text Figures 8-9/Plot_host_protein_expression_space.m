%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Analysis of the steady state interactions host-protein considering the mean value parameters
%   obtained from a set of runs of main_fitMu_ss_wildtype_mp5.m for Nr=57
%   as analysed in >Analysis_resultys_estim_Nr57_v5.m and Plot_results_average_Nr57_v5.m 
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

model_p.ku_A = 3.1120; %Dissotiation rate RBS-ribosome  (1/min)            
model_p.kb_A =  11.8005; %Assotiation rate RBS-ribosome  (1/min/molec)
model_p.Omega_A = 0.03; % Transcription rate   (1/min)           
model_p.lp_A = model_p.lp_nr; % Assume protein A with mean characteristics
model_p.dm_A = model_p.dm_nr;

proteinA_params =[model_p.ku_A, model_p.kb_A, model_p.Omega_A];

%%%%% EVALUATION FOR AVERAGE RBS AND VARYING PROMOTER:

varying_omega_A=[0.1*model_p.Omega_A, model_p.Omega_A,10*model_p.Omega_A,50*model_p.Omega_A,100*model_p.Omega_A,500*model_p.Omega_A,1000*model_p.Omega_A,1500*model_p.Omega_A,3000*model_p.Omega_A,5000*model_p.Omega_A];
for k=1:length(varying_omega_A)
    proteinA_params =[model_p.ku_A, model_p.kb_A,varying_omega_A(k)];
    
    [ke,KC0_r,KC0_nr, KC0_A,mu_estimated_profile,mp_estimated,mu_r_profile, JSum, JNr, J_A]  = eval_host_protein_ss(proteinA_params, f_substrate); 

    estimated_fraction_A =   J_A./JSum;
    estimated_fraction_R = JNr./JSum;
    estimated_fraction_P = 1-estimated_fraction_A-estimated_fraction_R;

    estimated_mass_A =   mp_estimated.*J_A./JSum;
    estimated_mass_productivity_A = estimated_mass_A.*mu_estimated_profile;
    estimated_specific_productivity_A = J_A./JSum.*mu_estimated_profile;
    
    r_estimated = mu_r_profile./mu_estimated_profile;
    rt_estimated = (1+JSum).*r_estimated/model_p.Phi_t;

    Struct_results_promoter.protein_chars{k}=proteinA_params;
    Struct_results_promoter.fraction_A{k}=estimated_fraction_A;
    Struct_results_promoter.fraction_R{k}=estimated_fraction_R;
    Struct_results_promoter.fraction_P{k}=estimated_fraction_P;
    Struct_results_promoter.mass_A{k}=estimated_mass_A;
    Struct_results_promoter.mass_productivity_A{k}=estimated_mass_productivity_A;
    Struct_results_promoter.specific_productivity_A{k}=estimated_specific_productivity_A;
    Struct_results_promoter.r{k} = r_estimated;
    Struct_results_promoter.rt{k} = rt_estimated;
    Struct_results_promoter.mu{k} = mu_estimated_profile;
end

figures_plot=true;
if figures_plot==true
ColorMap = jet(length(f_substrate))
Linewidth = 3;
Fontsize = 14;
    
     f1= figure(1)
         subplot(131)
         for m=1:length(f_substrate)
             line_fraction_A_si=[];
             line_fraction_P_si=[];
             line_fraction_R_si=[];
             for n=1:length(varying_omega_A)
               line_fraction_A_si=[line_fraction_A_si,Struct_results_promoter.fraction_A{n}(m)];
               line_fraction_P_si=[line_fraction_P_si,Struct_results_promoter.fraction_P{n}(m)];
               line_fraction_R_si=[line_fraction_R_si,Struct_results_promoter.fraction_R{n}(m)];
             end
             plot(log10(varying_omega_A),  line_fraction_A_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
              plot(log10(varying_omega_A),  line_fraction_R_si+line_fraction_P_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','Diamond', 'LineStyle','-','LineWidth',2);
         end
         grid on, xlabel('$\log_{10}{\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Phi_A (\mathrm{O}), \Phi_{P} + \Phi_{R} (\diamondsuit)$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         
         subplot(132)
         for m=1:length(f_substrate)
             line_mu_si=[];
             for n=1:length(varying_omega_A)
               line_mu_si=[line_mu_si,Struct_results_promoter.mu{n}(m)];
             end
             plot(log10(varying_omega_A), line_mu_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\mu  (\mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         
          subplot(133)
         for m=1:length(f_substrate)
             line_mass_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_mass_productivity_A_si=[ line_mass_productivity_A_si,Struct_results_promoter.mass_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  line_mass_productivity_A_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Pi_A  (g\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         
         colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = 'f(s_i)';
         hbc.Label.FontSize = 14;
          hbc.Title.FontSize = 12;
          
          %exportgraphics(f1,'./promoter_swap.png','Resolution',300)
          f2= figure(2)
           subplot(121)
         for m=1:length(f_substrate)
             line_mass_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_mass_productivity_A_si=[ line_mass_productivity_A_si,Struct_results_promoter.mass_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  line_mass_productivity_A_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\Pi_A  (g\cdot \mathrm{min}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         
          subplot(122)
         for m=1:length(f_substrate)
             line_specific_productivity_A_si=[];
             for n=1:length(varying_omega_A)
                line_specific_productivity_A_si=[ line_specific_productivity_A_si,Struct_results_promoter.specific_productivity_A{n}(m)];
             end
             plot(log10(varying_omega_A),  line_specific_productivity_A_si,'Color',ColorMap(m,:),'MarkerSize',10,...
                  'MarkerFaceColor',ColorMap(m,:), 'MarkerEdgeColor',[0 0 0],...
                  'Marker','o', 'LineStyle','-','LineWidth',2);
              hold on
         end
         grid on, xlabel('$\log_{10}{\omega_A}$ (mRNA/min)','FontSize',Fontsize,'Interpreter','latex')
         ylabel('$\pi_A  (g\cdot \mathrm{min}^{-1}\cdot \mathrm{cDW}^{-1})$','FontSize',Fontsize,'Interpreter','latex');
         hold off
         
         colormap(ColorMap);
         f_substrate_round=round(f_substrate, 2);
         hbc=colorbar('Direction','reverse');
         hbc.TickLabelsMode='manual';
         hbc.TickLabelInterpreter='latex';
         hbc.TickLabels=f_substrate_round;
         hbc.Title.String = 'f(s_i)';
         hbc.Label.FontSize = 14;
          hbc.Title.FontSize = 12;
         


end


