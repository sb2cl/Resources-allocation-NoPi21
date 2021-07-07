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
%global Best_Results_Nr57

num_runs=200;

Best_Results_Nr57=[];
for k=1:num_runs
 main_fitMu_ss_wildtype_mp_v6
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
Regular_average_results_sorted_choice_Nr57=mean(Best_Results_Nr57_sorted_choice,2)
Regular_std_results_sorted_choice_Nr57=std(Best_Results_Nr57_sorted_choice,1,2)


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

% Number of ribosomes (mature and inmature, active and inactive). molecs
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

model_p.Phi_h_b = 0.99*ones(5,1); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_h_t = 0.83*ones(5,1); %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_r_t = 0.3*ones(5,1);  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_nr_t = 0.55*ones(5,1);  %Initial iteration value (based on the results obtained from Analysis_results_estim_Nr57_v6.m)
model_p.Phi_m = 0.97; % Initial estimation of the fraction of mature to inmature ribosomes

model_p.model_mass=2; % Protein mass model  m =m0exp(βμ^gamma)  
model_p.mu_estimated = Bremer_exp_data.mu; %Initial estimations for the function eval_host_protein_ss_vScott.m
model_p.r_profile =  350*Bremer_exp_data.rt/Bremer_exp_data.rt(5); % molecs

%%%%% Weighted average values of estimated parameters

model_p.ku_r =  Weighted_average_results_Nr57(2,1); %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr = Weighted_average_results_Nr57(3,1); %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  Weighted_average_results_Nr57(4,1); %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  Weighted_average_results_Nr57(5,1); %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = Weighted_average_results_Nr57(6,1); % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr =Weighted_average_results_Nr57(7,1); % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.Phi_m = Weighted_average_results_Nr57(8,1); % Average ratio between mature available and total number of ribosomes 

%%%%% EVALUATION OF THE AVERAGE ESTIMATION:

evaluate_results=true;
if evaluate_results==true

[ke_r, ke_nr, KC0_r,KC0_nr,mu_estimated_profile,mh_estimated,mu_r_profile, JNr, JNnr]  = eval_Mu_ss_wildtype_mp_v6(Weighted_average_results_Nr57(2:8,1)'); 


    JWSum = model_p.WEm_r*JNr + model_p.WEm_nr *JNnr;
    JSum = JNr+ JNnr;
    Fraction_r=JNr./JSum; % Estimated fraction of ribosomal mass
    varphi_R = Bremer_exp_data.rt*model_p.ribosome_mass./Bremer_exp_data.mp; %Experimental fraction of ribosomal mass
    varphi_P = 1-varphi_R;
    estimated_varphi_R = Fraction_r;
    estimated_varphi_P = 1-estimated_varphi_R;

    r_estimated = mu_r_profile./mu_estimated_profile; %r_estimated = estimated_flux_mur/estimated_mu
    rt_estimated = (1+JWSum).*r_estimated/model_p.Phi_m; %Estimated total number of ribosomes: r_a = Phi_t*r_t -> r_t=r_a/Phi_t = r/(1+JWSum)/Phi_m

    r_mature_available = (1+JWSum).*r_estimated;
    r_active_estimated = model_p.Phi_h_t .*r_mature_available;  
    r_bound_at_RBS=r_mature_available-r_active_estimated;
    r_immature= rt_estimated-(r_mature_available+ r_estimated);
    
figures_plot=true;

if figures_plot==true
 
f1= figure(1)
subplot(131)
plot(Bremer_exp_data.mu, mu_estimated_profile,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineStyle','none'), 
xlabel('estimated $\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('experimental  $\mu (\mathrm{min}^{-1})$','Interpreter','latex'), grid on
hold on
fplot(@(x)x, [min(Bremer_exp_data.mu) max(Bremer_exp_data.mu)],'r-','Linewidth',2,'LineStyle','-')
ax = gca;
ax.FontSize = 14;
hold off
subplot(132)
plot(mu_estimated_profile, r_estimated,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$r$ (molec)','Interpreter','latex'), grid on
hold on
%plot(Bremer_exp_data.mu, 350*Bremer_exp_data.rt/Bremer_exp_data.rt(5),'r*','MarkerSize',10,'LineWidth',3)
%legend('Estimated free ribosomes', 'Initial estimation','Location','southeast')
ax = gca;
ax.FontSize = 14;
hold off
subplot(133)
plot(mu_estimated_profile, rt_estimated,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$r_t$ (molec)','Interpreter','latex'), grid on
hold on
plot(Bremer_exp_data.mu, Bremer_exp_data.rt ,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--')
legend('Estimated $r_t$', 'Experimental  $r_t$','Location','southeast','FontSize',14,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
hold off

%exportgraphics(f1,'./estimation_fractions.png','Resolution',300)


f2= figure(2)
subplot(121)
plot(mu_estimated_profile, ke_r,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
hold on
plot(mu_estimated_profile, ke_nr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
legend('$k_e^r$', '$k_e^{nr}$','Interpreter','latex','Location','northeast','FontSize',16)
xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$k_e$','Interpreter','latex'), grid on
ax = gca;
ax.FontSize = 14;
hold off
subplot(122)
plot(mu_estimated_profile,KC0_r,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--'), 
xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$K_{C0}^r$, $K_{C^0}^{nr}$','Interpreter','latex'), grid on
hold on
plot(mu_estimated_profile, KC0_nr,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--')
legend('$K_{C0}^r$', '$K_{C^0}^{nr}$','Interpreter','latex','Location','northeast','FontSize',16)
ax = gca;
ax.FontSize = 14;
hold off

%exportgraphics(f2,'../images_def/RBS_fractions_v6.png','Resolution',300)

 f3=figure(3);
plot(Bremer_exp_data.mu, varphi_R,'Marker','s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r'), grid on
hold on
plot(Bremer_exp_data.mu, varphi_P,'Marker','s','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b')
plot(mu_estimated_profile, estimated_varphi_R,'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','r') 
plot(mu_estimated_profile, estimated_varphi_P,'Marker','o','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',10,'LineWidth',2,'LineStyle','--','Color','b')

legend('Exp.$ \phi_R$','Exp. $\phi_P$','Extimated $\phi_R$','Estimated. $\phi_P$','Interpreter','latex','Location','southeast','FontSize',16)
 xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('$\phi_R$, $\phi_P$','Interpreter','latex')
 title(' Ribosomal ($\phi_R$) and non-ribosomal  ($\phi_P$) protein mass fractions','Interpreter','latex')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off
%exportgraphics(f3,'./mass_fractions.png','Resolution',300)


 f4=figure(4);
 subplot(121)
 R_s=[r_active_estimated'; r_bound_at_RBS' ;r_estimated';r_immature'];
area(mu_estimated_profile,R_s','LineWidth',2), grid on
legend('Active bound ribosomes', 'RBS bound ribosomes','Free ribosomes','Immature ribosomes', 'Location','northwest','FontSize',16)
 xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('number of ribosomes')
 title(' Estimated ribosomes')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off
 subplot(122)
 R_s=[r_active_estimated'./rt_estimated'; r_bound_at_RBS'./rt_estimated' ;r_estimated'./rt_estimated';r_immature'./rt_estimated'];
area(mu_estimated_profile,R_s','LineWidth',2), grid on
legend('Active bound ribosomes', 'RBS bound ribosomes','Free ribosomes','Immature ribosomes', 'Location','southeast','FontSize',16)
 xlabel('$\mu (\mathrm{min}^{-1})$','Interpreter','latex'), ylabel('fraction of ribosomes')
 title(' Estimated fractions of ribosomes')
 grid on
 ax = gca;
 ax.FontSize = 14;
hold off 
%exportgraphics(f4,'../images_def/r_ra_sums_v6.png','Resolution',300)

end
end

% Weighted_average_results_Nr57 =  mean of best 25 out of 200
%        OLD                           UPDATE                   STD UPDATE
%   345.1734                         346.0846                   2.2464
%   130.7432                        129.9034                    4.0728
%     3.1120                            3.0921                      0.1396
%     5.8138                            5.5745                      0.7809
%    11.8005                         12.8570                    1.4993
%     5.3606                          5.6509                        0.2914
%     0.0278                           0.02750                  0.0002552
%     0.7774                           0.9043                   0.005392


