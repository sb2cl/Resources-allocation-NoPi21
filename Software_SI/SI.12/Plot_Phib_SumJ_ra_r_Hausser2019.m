
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% Next we obtain the relationship between  Jk and the free ribosomes for each
% protein.
% We use the expressions:
%
%   Jk= lpk/nu*omegak*Kp/dmk/r
%
%   Emk = lpk/le*(1-(lpk/(lpk+le))^(lpk/le))
%
%  Kp: Experimental translation rates obtained from Hausser 2019
%  (protein/mRNA/time)
%  r: free ribosomes
%  lpk: protein length (aa)
%  le: 25 (aa) inverse of ribosome density
%  dmk= mRNA degradation rates (Hausser 2019 and Chen 2015)
%  omegaK: Experimental transcription rates obtained from Hausser 2019
%  (mRNA/time)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

[dades_temp,texte]=xlsread('./Dades/Hausser2019_coliRates/coliRates.xlsx','E. coli rates');
[x,y]=size(dades_temp);
[xt,yt]=size(texte);
Names_genes = texte(2:xt,1);

%Get columns of interest: translation rate, transcription rate, protein
%length, mRNA abundance, mRNA decay rate, protein abundance

Dades = [dades_temp(:,3),dades_temp(:,4),dades_temp(:,6),dades_temp(:,7) ...
  dades_temp(:,11),dades_temp(:,12)];

% For genes where there is no information on the mRNA degradation rate we
% consider the mean value in Che etal., Mol Syst Biol. 2015 Jan; 11(1):781:
 mean_dm= log(2)/2.5*60; % h^{-1}   
 for i=1:x
    if isnan(Dades(i,5))  
       Dades(i,5)=log10(mean_dm);
    end
 end

% We consider ribosomal proteins independently
Ribosomal_indices = 2058:2125;

% Get rid of genes at which some of the values bp,bm,lp is NaN
Bad_genes=any(isnan(Dades(:,1:3)),2);
Good_genes_all=find(1-Bad_genes); %Gives their indices 
Good_genes_ribosomal=Good_genes_all(ismember(Good_genes_all,Ribosomal_indices));
Good_genes_non_ribosomal = Good_genes_all(not(ismember(Good_genes_all,Ribosomal_indices')));

% Transcription, translation rates and protein length for all goodies:
Translation_rates_all=10.^Dades(Good_genes_all,1)/60;    %protein/(mRNA*min)
Length_proteins_all= 10.^Dades(Good_genes_all,3);         % amino acids (aa)
Transcription_rates_all=10.^Dades(Good_genes_all,2)/60;  %mRNA/min
Transcription_rates_all_aa= Transcription_rates_all.*Length_proteins_all; %aa/min

Mean_translation_rate_all = mean(Translation_rates_all)    %protein/(mRNA*min)
Mean_protein_length_all = mean(Length_proteins_all) %aa
Mean_transcription_rate_all = mean(Transcription_rates_all)    %mRNA/min
Mean_transcription_rate_all_aa = mean(Transcription_rates_all_aa)    %aa/min

% Transcription, translation rates and protein length for all ribosomal goodies:
Translation_rates_ribosomal=10.^Dades(Good_genes_ribosomal,1)/60;    %protein/(mRNA*min)
Length_proteins_ribosomal= 10.^Dades(Good_genes_ribosomal,3);         % amino acids (aa)
Transcription_rates_ribosomal=10.^Dades(Good_genes_ribosomal,2)/60;  %mRNA/min
Transcription_rates_ribosomal_aa= Transcription_rates_ribosomal.*Length_proteins_ribosomal; %aa/min
Names_ribosomes= Names_genes(Good_genes_ribosomal);
Total_length_proteins_ribosomal=sum(Length_proteins_ribosomal);

Mean_translation_rate_ribosomal = mean(Translation_rates_ribosomal)    %protein/(mRNA*min)
Mean_protein_length_ribosomal = mean(Length_proteins_ribosomal); %aa
Mean_transcription_rate_ribosomal = mean(Transcription_rates_ribosomal)    %mRNA/min
Mean_transcription_rate_ribosomal_aa = mean(Transcription_rates_ribosomal_aa)    %aa/min

% Transcription, translation rates and protein length for all NON ribosomal goodies:
Translation_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,1)/60;    %protein/(mRNA*min)
Length_proteins_non_ribosomal= 10.^Dades(Good_genes_non_ribosomal,3);         % amino acids (aa)
Transcription_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,2)/60;  %mRNA/min
Transcription_rates_non_ribosomal_aa= Transcription_rates_non_ribosomal.*Length_proteins_non_ribosomal; %aa/min

Mean_translation_rate_non_ribosomal = mean(Translation_rates_non_ribosomal)    %protein/(mRNA*min)
Mean_protein_length_non_ribosomal = mean(Length_proteins_non_ribosomal) %aa
Mean_transcription_rate_non_ribosomal = mean(Transcription_rates_non_ribosomal)    %mRNA/min
Mean_transcription_rate_non_ribosomal_aa = mean(Transcription_rates_non_ribosomal_aa)    %aa/min

% Mean mRNA degradation rate
mRNA_degradation_rates_all=10.^Dades(Good_genes_all,5)/60;  %1/min
mRNA_degradation_rates_ribosomal=10.^Dades(Good_genes_ribosomal,5)/60;  %1/min
mRNA_degradation_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,5)/60;  %1/min

Mean_dm_all = mean(mRNA_degradation_rates_all)
Mean_dm_non_ribosomal = mean(mRNA_degradation_rates_non_ribosomal)
Mean_dm_rate_ribosomal = mean(mRNA_degradation_rates_ribosomal)

%%%%%%% GENERAL CELL DATA:

mu_Hausser= log(2)/21.5; %0.032
maa = 182.6e-24; %average weight for an amino acid
% For the whole cell, Hausser 2019 estimates:
mc_Hausser = 180e-15;
nu=1260;
b=maa*nu/mc_Hausser; %b=5.3506e-07

td=20:5:100;
mu = log(2)./td;
Phib_t_ra = mu/b;
% Range of r assuming N(350,80) -4sigma to +2sigma
% R.J. Harvey, ﻿Fraction of ribosomes synthesizing protein as a function of
% specific growth rate. J. ﻿of Bacteriology, 1973. 
free_r = 30:1:510;

le=25;

rt_exp = [%72000
    45100
    26300
    13500
    6800];
Phi_t = 0.8;
ra_exp = Phi_t*rt_exp;
delta = (max(ra_exp)-min(ra_exp))/100;
Range_ra= min(ra_exp):delta:max(ra_exp);

%%%%%%% We get Jk for ALL, RIBOSOMAL and NON-RIBOSOMAL groups for Hausser

a_ratio_non_ribosomal=Transcription_rates_non_ribosomal_aa.*Translation_rates_non_ribosomal/nu./mRNA_degradation_rates_non_ribosomal;
Mean_a_ratio_non_ribosomal=mean(a_ratio_non_ribosomal)
Em_non_ribosomal = Length_proteins_non_ribosomal/le.*( 1-(Length_proteins_non_ribosomal./(Length_proteins_non_ribosomal+le)).^(Length_proteins_non_ribosomal/le) );
Mean_Em_non_ribosomal = mean(Em_non_ribosomal)   % We get 8.2431 which is extremely very close to  8.2415 obtained using the average lp=333
Em_weight_non_ribosomal = 1+ 1./Em_non_ribosomal;
Mean_Em_weight_non_ribosomal = mean(Em_weight_non_ribosomal) % We get 1.1929, slighly larger than 1+1/Mean_Em_non_ribosomal

a_ratio_ribosomal=Transcription_rates_ribosomal_aa.*Translation_rates_ribosomal/nu./mRNA_degradation_rates_ribosomal;
Mean_a_ratio_ribosomal=mean(a_ratio_ribosomal)
Em_ribosomal = Length_proteins_ribosomal/le.*( 1-(Length_proteins_ribosomal./(Length_proteins_ribosomal+le)).^(Length_proteins_ribosomal/le) );
Mean_Em_ribosomal = mean(Em_ribosomal )
Em_weight_ribosomal = 1+ 1./Em_ribosomal;
Mean_Em_weight_ribosomal = mean(Em_weight_ribosomal) % We get 1.3566

a_ratio_all=Transcription_rates_all_aa.*Translation_rates_all/nu./mRNA_degradation_rates_all;
Mean_a_ratio_all=mean(a_ratio_all)
Em_all = Length_proteins_all/le.*( 1-(Length_proteins_all./(Length_proteins_all+le)).^(Length_proteins_all/le) );
Mean_Em_all = mean(Em_all )
Em_weight_all = 1+ 1./Em_all;
Mean_Em_weight_all = mean(Em_weight_all) % We get  1.1960

% The sum of "Jks" and weighted Jks
Sum_Jk_non_ribosomal=0;
Sum_Jk_weighted_non_ribosomal=0;
for k=1:size(Translation_rates_non_ribosomal,1)
    Sum_Jk_non_ribosomal=Sum_Jk_non_ribosomal+a_ratio_non_ribosomal(k); 
    Sum_Jk_weighted_non_ribosomal=Sum_Jk_weighted_non_ribosomal+a_ratio_non_ribosomal(k)*Em_weight_non_ribosomal(k); 
end
Mean_Jk_non_ribosomal=Sum_Jk_non_ribosomal/size(Translation_rates_non_ribosomal,1)
Mean_Jk_weighted_non_ribosomal=Sum_Jk_weighted_non_ribosomal/size(Translation_rates_non_ribosomal,1)
Sum_Jk_non_ribosomal
Sum_Jk_weighted_non_ribosomal

Sum_Jk_ribosomal=0;
Sum_Jk_weighted_ribosomal=0;
for k=1:size(Translation_rates_ribosomal,1)
    Sum_Jk_ribosomal=Sum_Jk_ribosomal+a_ratio_ribosomal(k); 
    Sum_Jk_weighted_ribosomal=Sum_Jk_weighted_ribosomal+a_ratio_ribosomal(k)*Em_weight_ribosomal(k); 
end
Mean_Jk_ribosomal=Sum_Jk_ribosomal/size(Translation_rates_ribosomal,1)
Mean_Jk_weighted_ribosomal=Sum_Jk_weighted_ribosomal/size(Translation_rates_ribosomal,1)
Sum_Jk_ribosomal
Sum_Jk_weighted_ribosomal

Sum_Jk_all=0;
Sum_Jk_weighted_all=0;
for k=1:size(Translation_rates_all,1)
    Sum_Jk_all=Sum_Jk_all+a_ratio_all(k);  
    Sum_Jk_weighted_all=Sum_Jk_weighted_all+a_ratio_all(k)*Em_weight_all(k); 
end
Mean_Jk_all=Sum_Jk_all/size(Translation_rates_all,1)
Mean_Jk_weighted_all=Sum_Jk_weighted_all/size(Translation_rates_all,1)
Sum_Jk_all
mu_Hausser/(b*Sum_Jk_all) % Must give 1. From Hausser data we get 1.151331124685091
Phib_t_ra
Sum_Jk_non_ribosomal+ Sum_Jk_ribosomal
Sum_Jk_weighted_all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We obtain Phi_b_t, Phi_b as a function of r

ra_Hausser = free_r +  Sum_Jk_weighted_all;

Phi_p_t_Hausser = Sum_Jk_non_ribosomal./(free_r +Sum_Jk_weighted_all);

Phi_r_t_Hausser = Sum_Jk_ribosomal./(free_r +Sum_Jk_weighted_all);

Phi_t_b_Hausser = Phi_p_t_Hausser  + Phi_r_t_Hausser;

Phi_p_Hausser = Sum_Jk_weighted_non_ribosomal./(free_r +Sum_Jk_weighted_all);

Phi_r_Hausser = Sum_Jk_weighted_ribosomal./(free_r +Sum_Jk_weighted_all);

Phi_b_Hausser = Sum_Jk_weighted_all./(free_r +Sum_Jk_weighted_all);

Phib_t_ra_Hausser = mu_Hausser/b; 

[R,RA] = meshgrid(free_r,Range_ra);
PHIB = (RA-R)./RA;
SUMJ_weighted = PHIB./(1-PHIB);

Fontsize=14;
f1=figure(1)
subplot(121)
 contour(R,RA,PHIB,'ShowText','on','LineWidth',3);
 hold on
 plot(free_r,ra_Hausser,'k','LineWidth',5)
 xlabel('r'), ylabel('r_a')
 title('\Phi_b','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 14;
% ax.YLim = [0.01 0.16];

 subplot(122)
contour(R,RA,SUMJ_weighted,'ShowText','on','LineWidth',3)
 xlabel('r'), ylabel('r_a')
 title('$\sum_{k=\left\{p,r\right\}}[1+\frac{1}{E_{mk}}]J_k(\mu, r)$','Interpreter','latex','FontSize',Fontsize)
 grid on
 hold on
 plot(free_r,ra_Hausser,'k','LineWidth',5)
 ax = gca;
 ax.FontSize = 14;
% ax.YLim = [0.01 0.16];
%print -dpng plot_Sum_J_r_ra.png
%f1.WindowState = 'maximized';
%exportgraphics(f1,'../images_def/plot_Phib_Sum_J_r_ra.png','Resolution',300)

Fontsize=16;

f2=figure(2)
 plot(free_r,Phi_b_Hausser,'LineWidth',3,'Color','k')
 hold on
 plot(free_r,Phi_t_b_Hausser,'LineWidth',3,'Color','b')
 legend('$\Phi_b$',' $\Phi^t_b$','Interpreter','latex','FontSize',Fontsize) 
 xlabel('r','Interpreter','latex','FontSize',Fontsize), ylabel('$\Phi_b, \Phi^t_b$','Interpreter','latex','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 16;
% ax.YLim = [0.01 0.16];
%print -dpng ../Texte/images/f_Sum_J.png
%f2.WindowState = 'maximized';
%exportgraphics(f1,'../images_def/plot_Phib_Phibt_r.png','Resolution',300)


 
 
