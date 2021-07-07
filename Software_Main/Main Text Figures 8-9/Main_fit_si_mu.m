

Bremer_exp_data.f_si = [
    12
    16
    18
    20
    21]*60/1260; 

Bremer_exp_data.mu = [
    log(2)/100
    log(2)/60
    log(2)/40
    log(2)/30
    log(2)/24]; %1/min


global Bremer_exp_data;

problem2.f='Cost_fit_si_mu';    %script with the cost function to optimize
problem2.x_L=[0,10,10,0,0]; % minimum expected values for b0,b1,a1
problem2.x_U=[0.001,20,35,2,2]; %  maximum expected values
opts2.maxeval=7500; 
Results2=MEIGO(problem2,opts2,'ESS'); 
% Best solution value		0.000630984
% Decision vector
% 	4.45687e-05
% 	19.9843
% 	34.9866
% 	0.679116
% 	1.06627
 f1=figure(1);
  plot(Bremer_exp_data.mu,Bremer_exp_data.f_si,'b*','MarkerSize',10,'LineWidth',3)
  hold on
  fplot(@(x)(Results2.xbest(1,1) + Results2.xbest(1,2)*x.^Results2.xbest(1,4))./(1+Results2.xbest(1,3)*x.^Results2.xbest(1,5)),[0,0.03],'r-','LineWidth',3)
  xlabel('\mu (min^{-1})'), ylabel('f(s_i)')
  grid on
 ax = gca;
 ax.FontSize = 13;


td_points=24:6:100;
mu_points= log(2)./td_points;

f_si_points= (Results2.xbest(1,1) + Results2.xbest(1,2)*mu_points.^Results2.xbest(1,4))./(1+Results2.xbest(1,3)*mu_points.^Results2.xbest(1,5));

plot(mu_points,f_si_points,'g*','MarkerSize',10,'LineWidth',3)
hold off
