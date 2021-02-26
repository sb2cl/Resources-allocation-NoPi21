%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Evaluation of the resources recruitment strengths (J_N_r, J_N_nr, J_A) at
% steady state  for the wild-type model 
%  
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ke,KC0_r,KC0_nr, KC0_A, mu_estimated_profile,mp_estimated,mu_r_profile, JSum, JNr, J_A] = eval_host_protein_ss(proteinA_params, f_substrate) 

global model_p;

ku_A = proteinA_params(1,1);
kb_A = proteinA_params(1,2);
Omega_A = proteinA_params(1,3);

nu_t_si = f_substrate.*model_p.nu_max; % vector with nu_t(s_i) 
ke= nu_t_si/model_p.l_e;
KC0_r =model_p.kb_r./(model_p.ku_r+ ke);
KC0_nr =model_p.kb_nr./(model_p.ku_nr+ ke);
KC0_A =kb_A./(ku_A+ ke);

while true 
    
if model_p.model_mass==1
 m0=77.3748e-15; %(g)  
  beta=61.7813; % (min)
 %mp_estimated = m0*exp(beta*mu_exp); % Notice that uses the experimental mu, not the estimated one
 mp_estimated = m0*exp(beta*model_p.mu_estimated); 
else
    m0=1.29181e-14;
    beta=	14.1089;
    gamma= 0.389004;
   %mp_estimated = m0*exp(beta*mu_exp.^gamma);
    mp_estimated = m0*exp(beta*model_p.mu_estimated.^gamma);
end

mu_r =  model_p.m_aa.*nu_t_si.*mp_estimated.*(1-model_p.Phi_b)./model_p.Phi_b.*(model_p.Phi_t.*model_p.Phi_r/model_p.ribosome_mass).^2;

%while true 
     
    J_N_r_estimated = 0.62*model_p.N_r*model_p.lp_r/model_p.l_e*model_p.Omega_r./(model_p.dm_r./KC0_r + mu_r);
    J_N_nr_estimated =  0.62*model_p.N_nr*model_p.lp_nr/model_p.l_e*model_p.Omega_nr./(model_p.dm_nr./KC0_nr + mu_r);
    J_A_estimated =  0.62*model_p.lp_A/model_p.l_e*Omega_A./(model_p.dm_A./KC0_A+ mu_r);

     Jsum=J_N_r_estimated + J_N_nr_estimated+J_A_estimated;
     model_p.Phi_r = J_N_r_estimated./(1+Jsum);
     model_p.Phi_b = Jsum./(1+Jsum);
     
     mu_r_new =  model_p.m_aa.*nu_t_si.*mp_estimated.*(1-model_p.Phi_b)./model_p.Phi_b.*(model_p.Phi_t.*model_p.Phi_r/model_p.ribosome_mass).^2;
       
    if sum(abs(mu_r_new - mu_r)) < 1e-3
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_si.*model_p.Phi_t.*model_p.Phi_r;
        mu_r_profile =  mu_r;
        model_p.mu_estimated = mu_estimated_profile;
        break;
    else
        mu_r = mu_r_new;
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_si.*model_p.Phi_t.*model_p.Phi_r;
        mu_r_profile =  mu_r;
        model_p.mu_estimated = mu_estimated_profile;
    end
end

JSum=model_p.Phi_b./(1-model_p.Phi_b);
JNr=model_p.Phi_r.*(1+JSum);
J_A=J_A_estimated;

end

  



   
     
   
 
 
 