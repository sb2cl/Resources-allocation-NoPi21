

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Estimation of the steady state values for the wild-type model 
%
% In this script we use the phenomenological relationship mp(mu) obtained
%  from thge experimental data:
%   m =m0exp(βμ)    m0=77.3748e−15(g)    β=61.7813(min)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Global variables:
global model_p;
global Bremer_exp_data;

% We fix the initial protein mass. We use the values in Bremer (cite)
% corresponding to the growth rates
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

model_p.l_e = 25; 	%Ribosome occupancy length (aa)
model_p.lp_r = 195; %'Mean length of ribosomal proteins (aa)     
model_p.lp_nr = 333; %'Mean length of non-ribosomal proteins (aa)     
model_p.nu_max = 1260; % Maximum effective translation rate per ribosome (aa/min)
model_p.dm_r  =  0.16; %'Mean degradation rate of ribosomal mRNA (1/min)
model_p.dm_nr  =  0.2; %'Mean degradation rate of non-ribosomal mRNA (1/min)
model_p.N_r = 57;%56;% 57; %Number of protein types building up a ribosome (molec)			
model_p.N_nr = 1735; %Number of non ribosomal protein types being expressed at one given time (molec)
model_p.WEm_r = 1.2891; % Average weight (1+1/Emr) for ribosomal protein-coding genes
model_p.WEm_nr =  1.1575; % Average weight (1+1/Emp) for non-ribosomal protein-coding genes
model_p.m_aa = 182.6e-24; % Average amino acid mass (g/aa)
% Dalton unit
Da = 1.6605e-24; % g/Da
%model_p.ribosome_mass = 2.7e6*Da; % Ribosome weight (g)
%model_p.ribosome_mass = model_p.m_aa*model_p.N_r*model_p.lp_r; 
model_p.ribosome_mass = 2.29*model_p.m_aa*model_p.N_r*model_p.lp_r; % mp_estimated includes the RNA component of ribosomes.  %CHECK

% Initial estimation of model parameters to be obtained:

model_p.ku_r =  132.68; %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr =  6; %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  7.82; %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  15; %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = 2.4; % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = 0.05; % Average transcription rate for non-ribosomal proteins  (1/min)        
model_p.r_profile =  350*Bremer_exp_data.rt/Bremer_exp_data.rt(5); % molecs
model_p.Phi_b = 0.99*ones(5,1); %Initial estimation (based on the results obtained from Hausser's data for r=1)
model_p.Phi_b_t = 0.8*ones(5,1); %Initial estimation (based on the results obtained from Hausser's data for r=1)
model_p.Phi_r_t = 0.3*ones(5,1);  %Initial estimation (based on the results obtained from Hausser's data for r=1)
model_p.Phi_p_t = 0.5*ones(5,1);  %Initial estimation (based on the results obtained from Hausser's data for r=1)
model_p.Phi_t = Bremer_exp_data.Phi_t ; 
model_p.mu_estimated = Bremer_exp_data.mu;
%model_p.model_mass=1; % Protein mass model  m =m0exp(βμ)    m0=77.3748e−15(g)    β=61.7813(min)
model_p.model_mass=2; % Protein mass model  m =m0exp(βμ^gamma)  
            
%========================= PROBLEM SPECIFICATIONS ===========================
problem.f='CostF_Mu_ss_wildtype_mp_v6';    %script with the cost function to optimize

problem.x_L=[100,3,3,3,4,0.02,0.77]; % minimum expected values for ku_r, ku_nr, kb_r, kb_nr,omega_r, omega_nr, Phi_t
problem.x_U=[135,10,15,15,6,0.04,0.999]; %  maximum expected values, Phi_t

opts.maxeval=7500; 
%opts.maxtime=1e7;  
%opts.local.solver='dhc';
%opts.local.iterprint=0;
%opts.local.n1=2;
%opts.local.n2=3;
%========================= END OF PROBLEM SPECIFICATIONS =====================

Results=MEIGO(problem,opts,'ESS'); 
vpa(Results.xbest')

model_p.ku_r =  Results.xbest(1,1); %Dissotiation rate RBS-ribosome for ribosomal mRNA (1/min)
model_p.ku_nr =  Results.xbest(1,2); %Dissotiation rate RBS-ribosome for non-ribosomal mRNA (1/min)            
model_p.kb_r =  Results.xbest(1,3); %Assotiation rate RBS-ribosome for ribosomal mRNA (1/min/molec)
model_p.kb_nr =  Results.xbest(1,4); %Assotiation rate RBS-ribosome for non-ribosomal mRNA (1/min/molec)                	
model_p.Omega_r = Results.xbest(1,5); % Average transcription rate for ribosomal proteins  (1/min)           
model_p.Omega_nr = Results.xbest(1,6); % Average transcription rate for non-ribosomal proteins  (1/min)      
model_p.Phi_t = Results.xbest(1,7); % Fraction of available w.r.t total ribosomes     

%%%%%%% PLOTTING RESULTS  %%%%%%%%%%%

evaluate_results=true;
if evaluate_results==true

[ke,KC0_r,KC0_nr, mu_estimated_profile,mp_estimated,mu_r_profile] = eval_Mu_ss_wildtype_mp_v6(Results.xbest); 

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

Results.xbest
figures_plot=false;

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

%Nr=57: 131.0464    3.1564    6.9172   12.0663    4.9323    0.0297    0.8002




