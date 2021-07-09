%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Evaluation of the resources recruitment strengths (J_N_r, J_N_nr, J_A) at
% steady state  for the average wild-type model and a protein of interest A
%  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ke_A,KC0_r,KC0_nr, KC0_A,mu_estimated_profile,estimated_mass_h,mu_r_profile, JNr, JNnr, JNA, JWSum] = eval_host_protein_ss_v7(proteinA_params, f_substrate) 

global model_p;

ku_A = proteinA_params(1,1);
kb_A = proteinA_params(1,2);
Omega_A = proteinA_params(1,3);

nu_t_si = f_substrate.*model_p.nu_max; % row vector with nu_t(s_i) 
ke_A= nu_t_si/model_p.le_A;
ke_r= nu_t_si/model_p.le_r;
ke_nr= nu_t_si/model_p.le_nr;
KC0_r =model_p.kb_r./(model_p.ku_r+ ke_r);
KC0_nr =model_p.kb_nr./(model_p.ku_nr+ ke_nr);
KC0_A =kb_A./(ku_A+ ke_A);

while true 
    
if model_p.model_mass==1
 m0=77.3748e-15; %(g)  
  beta=61.7813; % (min)
 %mp_estimated = m0*exp(beta*mu_exp); % Notice that uses the experimental mu, not the estimated one
 estimated_mass_h = m0*exp(beta*model_p.mu_estimated); 
else
    m0=1.29181e-14;
    beta=	14.1089;
    gamma= 0.389004;
   %mp_estimated = m0*exp(beta*mu_exp.^gamma);
    estimated_mass_h = m0*exp(beta*model_p.mu_estimated.^gamma); % this is the mass of the native host
end

% We first estimate the flux mu*r of free resources:
     mu_r =  model_p.m_aa.*nu_t_si.*estimated_mass_h.*(1-model_p.Phi_s_b)./model_p.Phi_h_t.*(model_p.Phi_m .*model_p.Phi_r_t/model_p.ribosome_mass).^2;

%Now we estimate the sums Nr*Jr, Np*Jp  and J_A

    JNr_estimated =  model_p.N_r*model_p.Em_r*model_p.Omega_r./(model_p.dm_r./KC0_r + mu_r);
    JNnr_estimated =  model_p.N_nr*model_p.Em_nr*model_p.Omega_nr./(model_p.dm_nr./KC0_nr + mu_r);
    JNA_estimated =  model_p.Em_A*Omega_A./(model_p.dm_A./KC0_A+ mu_r);

    JWSum = model_p.WEm_r*JNr_estimated + model_p.WEm_nr *JNnr_estimated + model_p.WEm_A*JNA_estimated;
    
     model_p.Phi_s_b = JWSum./(1+JWSum);
     model_p.Phi_r_t = JNr_estimated./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for ribosomal protein-coding genes
     model_p.Phi_nr_t = JNnr_estimated./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for non-ribosomal endogenous  protein-coding genes
     model_p.Phi_h_t =  model_p.Phi_r_t  + model_p.Phi_nr_t;

    % We update the estimate of the flux mu*r of free resources:
     mu_r_new =  model_p.m_aa.*nu_t_si.*estimated_mass_h.*(1-model_p.Phi_s_b)./model_p.Phi_h_t.*(model_p.Phi_m .*model_p.Phi_r_t/model_p.ribosome_mass).^2;

    if sum(abs(mu_r_new - mu_r)) < 1e-3
         % We update the estimation of the cell growth rate:
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_si.*model_p.Phi_m.*model_p.Phi_r_t;
        mu_r_profile =  mu_r;
        model_p.mu_estimated = mu_estimated_profile;
        break;
    else
        mu_r = mu_r_new;
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_si.*model_p.Phi_m .*model_p.Phi_r_t;
        mu_r_profile =  mu_r;
        model_p.mu_estimated = mu_estimated_profile;
    end
end %while

% Update values
JNr= model_p.N_r*model_p.Em_r*model_p.Omega_r./(model_p.dm_r./KC0_r + mu_r_profile);
JNnr= model_p.N_nr*model_p.Em_nr*model_p.Omega_nr./(model_p.dm_nr./KC0_nr + mu_r_profile);
JNA =  model_p.Em_A*Omega_A./(model_p.dm_A./KC0_A+ mu_r_profile);

JWSum = model_p.WEm_r*JNr + model_p.WEm_nr *JNnr + model_p.WEm_A*JNA;
model_p.Phi_s_b = JWSum./(1+JWSum);
model_p.Phi_r_t = JNr./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for ribosomal protein-coding genes
model_p.Phi_nr_t = JNnr./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for non-ribosomal endogenous  protein-coding genes
model_p.Phi_h_t =  model_p.Phi_r_t  + model_p.Phi_nr_t;

if model_p.model_mass==1
     estimated_mass_h = m0*exp(beta*model_p.mu_estimated); 
else
    estimated_mass_h = m0*exp(beta*model_p.mu_estimated.^gamma); % this is the mass of the native host
end

end

  



   
     
   
 
 
 