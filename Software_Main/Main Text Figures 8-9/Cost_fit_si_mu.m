function J=Cost_fit_si_mu(x)
global Bremer_exp_data
mu= Bremer_exp_data.mu;
f_si= Bremer_exp_data.f_si;

  b_0 =  x(1,1);
  b_1 =  x(1,2);
  a_1 = x(1,3);
  c_1 = x(1,4);
  c_2 = x(1,5);
  prediction_f_si = (b_0 + b_1*mu.^c_1)./(1+ a_1*mu.^c_2);
  e= f_si-prediction_f_si ;
  J = e'*e;
end