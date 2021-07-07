function J=mse_wildtype(x)
global data_log10mur
global data_log10r
datax = data_log10mur;
datay = data_log10r;
  b_mse =  x(1,1);
  k_mse =  x(1,2);
  prediction = b_mse + k_mse*datax;
  e= datay-prediction;
  J = e'*e;
end