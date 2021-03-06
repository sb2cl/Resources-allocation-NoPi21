%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Evaluation of the resources recruitment strengths (J_N_r, J_N_nr) at
% steady state  for the wild-type model 
%  
%  The function eval_J_ss_wildtype has input arguments:
%       x = set of parameters kb_r, ku_r, kb_nr, ku_nr, omega_r, omega_nr
%       being optimized
%       r = profile r(mu) from the previous optimization iteration
%
%   In this script we use the phenomenological relationship mp(mu) obtained
%  from thge experimental data:
%   m =m0*exp(βμ)    m0=77.3748e−15(g)    β=61.7813(min)
%
%  WARNING: notice we use the subindexes p and nr interchangeably  to denote
%  non-ribosomal proteins
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ke,KC0_r,KC0_nr, mu_estimated_profile,mp_estimated,mu_r_profile] = eval_Mu_ss_wildtype_mp_v6(estimated_params) 

global model_p;
global Bremer_exp_data;

ku_r = estimated_params(1,1);
ku_nr = estimated_params(1,2);
kb_r = estimated_params(1,3);
kb_nr = estimated_params(1,4);
Omega_r = estimated_params(1,5);
Omega_nr = estimated_params(1,6);
Phi_t_mod = estimated_params(1,7);

nu_t_exp = Bremer_exp_data.nu_t; % vector with nu_t for each mu
mu_exp =  Bremer_exp_data.mu;
Phi_t_exp = Bremer_exp_data.Phi_t;
ke= nu_t_exp/model_p.l_e;
KC0_r =kb_r./(ku_r+ ke);
KC0_nr =kb_nr./(ku_nr+ ke);

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

% We first estimate the flux mu*r of free resources:
mu_r =  model_p.m_aa.*nu_t_exp.*mp_estimated.*(1-model_p.Phi_b)./model_p.Phi_b_t.*(Phi_t_mod.*model_p.Phi_r_t/model_p.ribosome_mass).^2;

%Now we estimate the sums Nr*Jr and Np*Jp 
    %a_nr = 1-(model_p.lp_nr/(model_p.lp_nr+model_p.l_e))^(model_p.lp_nr/model_p.l_e); %Use instead of factor 0.62 if more precision is wanted. Gives 0.6187
    %a_r = 1-(model_p.lp_r/(model_p.lp_r+model_p.l_e))^(model_p.lp_r/model_p.l_e); %Use instead of factor 0.62 if more precision is wanted. Gives 0.6097

     
     J_N_r_estimated = 0.62*model_p.N_r*model_p.lp_r/model_p.l_e*Omega_r./(model_p.dm_r./KC0_r + mu_r);
    J_N_nr_estimated =  0.62*model_p.N_nr*model_p.lp_nr/model_p.l_e*Omega_nr./(model_p.dm_nr./KC0_nr + mu_r);

     JSum=J_N_r_estimated + J_N_nr_estimated;
     JWSum = model_p.WEm_r*J_N_r_estimated + model_p.WEm_nr *J_N_nr_estimated;
     
     model_p.Phi_b = JWSum./(1+JWSum);
     model_p.Phi_r_t = J_N_r_estimated./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for ribosomal protein-coding genes
     model_p.Phi_p_t = J_N_nr_estimated./(1+JWSum);  %Fraction of actively translating ribosomes w.r.t. mature available ones for non-ribosomal protein-coding genes
     model_p.Phi_b_t =  model_p.Phi_r_t  + model_p.Phi_p_t;

 % We update the estimate of the flux mu*r of free resources:
     mu_r_new =  model_p.m_aa.*nu_t_exp.*mp_estimated.*(1-model_p.Phi_b)./model_p.Phi_b_t.*(Phi_t_mod.*model_p.Phi_r_t/model_p.ribosome_mass).^2;
       
    if sum(abs(mu_r_new - mu_r)) < 1e-3
        % We update the estimation of the cell growth rate:
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_exp.*Phi_t_mod.*model_p.Phi_r_t;
        model_p.mu_estimated = mu_estimated_profile;
        break;
    else
         mu_r = mu_r_new;
        mu_estimated_profile = model_p.m_aa./model_p.ribosome_mass.*nu_t_exp.*Phi_t_mod.*model_p.Phi_r_t;
        mu_r_profile =  mu_r;
        model_p.mu_estimated = mu_estimated_profile;
    end
end

end

  



   
     
   
 
 
 