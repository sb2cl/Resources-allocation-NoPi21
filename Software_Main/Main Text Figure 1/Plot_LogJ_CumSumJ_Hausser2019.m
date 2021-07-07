
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      
% Next we obtain the MEAN relationship between  Jk and the free ribosomes for each
% protein.
% We use the expression:
%
%   Jk= lpk/nu*omegak*Kp/dmk/r
%
%  Kp: Experimental translation rates obtained from Hausser 2019
%  (protein/mRNA/time)
%  r: free ribosomes
%  lpk: protein length (aa)
%  dmk= mRNA degradation rates obtained from Hausser 2019 (1/time)
%  omegaK: Experimental transcription rates obtained from Hausser 2019
%  (mRNA/time)
%
% We consider the logarithmic relationship:
%
%  log10(Jk) = log10(lpk/nu*omegak*Kp/dmk) - log10(r)
%
% and then we get the mean for the different groups of proteins
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

Avg99_non_ribosomal= Sum_Jk_non_ribosomal/1735
Avg95_non_ribosomal= Sum_Jk_non_ribosomal/875
Avg99_weighted_non_ribosomal= Sum_Jk_weighted_non_ribosomal/1735
Avg95_weighted_non_ribosomal= Sum_Jk_weighted_non_ribosomal/875

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

Avg99_ribosomal= Sum_Jk_ribosomal/57
Avg95_ribosomal= Sum_Jk_ribosomal/49
Avg99_weighted_ribosomal= Sum_Jk_weighted_ribosomal/57
Avg95_weighted_ribosomal= Sum_Jk_weighted_ribosomal/49

%%% OBTAINING AVERAGE VALUES:

Mean_J_non_ribosomal_99 = Sum_Jk_non_ribosomal/1735
Mean_J_ribosomal_99 = Sum_Jk_ribosomal/57

Mean_Em_weight_ribosomal_99 = Sum_Jk_weighted_ribosomal/Sum_Jk_ribosomal . % 1.289082535135112
Mean_Em_weight_non_ribosomal_99 = Sum_Jk_weighted_non_ribosomal/Sum_Jk_non_ribosomal . % 1.157480361422457


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


%%%%%%% We ORDER the Jks for ALL, RIBOSOMAL and NON-RIBOSOMAL groups 
% [B,I]=sort(A)
% B lists the sorted dates and I contains the corresponding indices of A.
[a_ratio_non_ribosomal_sorted,a_ratio_non_ribosomal_sorted_index]=sort(a_ratio_non_ribosomal,'descend');
cumsum_a_ratio_non_ribosomal_sorted = cumsum(a_ratio_non_ribosomal_sorted);
[a_ratio_ribosomal_sorted,a_ratio_ribosomal_sorted_index]=sort(a_ratio_ribosomal,'descend');
cumsum_a_ratio_ribosomal_sorted = cumsum(a_ratio_ribosomal_sorted);
[a_ratio_all_sorted,a_ratio_all_sorted_index]=sort(a_ratio_all,'descend');
cumsum_a_ratio_all_sorted = cumsum(a_ratio_all_sorted);
who_is_ribosomal_index=ismember(Good_genes_all,Ribosomal_indices');
 
 f1 = figure(1);
subplot(221)
plot(a_ratio_non_ribosomal_sorted,'b.','MarkerSize',10)
 xlabel('Non-ribosomal proteins'), ylabel('J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
subplot(222)
plot(a_ratio_ribosomal_sorted,'b.','MarkerSize',10)
 xlabel('Ribosomal proteins'), ylabel('J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
%  subplot(233)
% plot(a_ratio_all_sorted,'b.','MarkerSize',10)
%  xlabel('All proteins'), ylabel('J_k(r=1) ')
%  grid on
%  ax = gca;
%  ax.FontSize = 13;
 subplot(223)
plot(a_ratio_non_ribosomal_sorted./Length_proteins_non_ribosomal(a_ratio_non_ribosomal_sorted_index),'b.','MarkerSize',10)
 xlabel('Non-ribosomal proteins'), ylabel('J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
subplot(224)
plot(a_ratio_ribosomal_sorted./Length_proteins_ribosomal(a_ratio_ribosomal_sorted_index),'b.','MarkerSize',10)
 xlabel('Ribosomal proteins'), ylabel('J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
%  subplot(236)
% plot(a_ratio_all_sorted./Length_proteins_all(a_ratio_all_sorted_index),'b.','MarkerSize',10)
%  xlabel('All proteins'), ylabel('J_k/l_k(r=1) ')
%  grid on
%  ax = gca;
%  ax.FontSize = 13;
% print -dpng ../Texte/images/J_sorted.png

%  f1.WindowState = 'maximized';
% exportgraphics(f1,'../Texte/images/J_sorted.png','Resolution',300)

f2 = figure(2);
subplot(121)
plot(cumsum_a_ratio_non_ribosomal_sorted,'b.','MarkerSize',15)
hold on
fplot(@(x)0.95*cumsum_a_ratio_non_ribosomal_sorted(end),[1 length(cumsum_a_ratio_non_ribosomal_sorted)],'r-')
fplot(@(x)0.99*cumsum_a_ratio_non_ribosomal_sorted(end),[1 length(cumsum_a_ratio_non_ribosomal_sorted)],'g-')
%fplot(0.99*cumsum_a_ratio_non_ribosomal_sorted(end)*ones(length(cumsum_a_ratio_non_ribosomal_sorted)),'g-')
 xlabel('Non-ribosomal proteins'), ylabel('CumSum J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
 hold off
subplot(122)
plot(cumsum_a_ratio_ribosomal_sorted,'b.','MarkerSize',15)
hold on
fplot(@(x)0.95.*cumsum_a_ratio_ribosomal_sorted(end),[1 length(cumsum_a_ratio_ribosomal_sorted)],'r-')
fplot(@(x)0.99.*cumsum_a_ratio_ribosomal_sorted(end),[1 length(cumsum_a_ratio_ribosomal_sorted)],'g-')
%plot(0.95*cumsum_a_ratio_ribosomal_sorted(end)*ones(length(cumsum_a_ratio_ribosomal_sorted)),'r-')
%plot(0.99*cumsum_a_ratio_ribosomal_sorted(end)*ones(length(cumsum_a_ratio_ribosomal_sorted)),'g-')
 xlabel('Ribosomal proteins'), ylabel('CumSum J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
 hold off
%  subplot(133)
% plot(cumsum_a_ratio_all_sorted,'b.','MarkerSize',10)
% hold on
% fplot(@(x)0.95.*cumsum_a_ratio_all_sorted(end),[1 length(cumsum_a_ratio_all_sorted)],'r-')
% fplot(@(x)0.99.*cumsum_a_ratio_all_sorted(end),[1 length(cumsum_a_ratio_all_sorted)],'g-')
% %plot(0.95*cumsum_a_ratio_all_sorted(end)*ones(length(cumsum_a_ratio_all_sorted)),'r-')
% %plot(0.99*cumsum_a_ratio_all_sorted(end)*ones(length(cumsum_a_ratio_all_sorted)),'g-')
%  xlabel('All proteins'), ylabel('CumSum J_k(r=1) ')
%  grid on
%  ax = gca;
 %ax.FontSize = 13;
 %print -dpng ../Texte/images/CumSum_J_sorted.png
 
%  f2.WindowState = 'maximized';
%  exportgraphics(f2,'../Texte/images/CumSum_J_sorted.png','Resolution',300)
 
f3 = figure(3);
subplot(121)
plot(Length_proteins_all,'b.','MarkerSize',10)
 xlabel('Non-ribosomal proteins'), ylabel('J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
 subplot(122)
plot(Length_proteins_all(a_ratio_all_sorted_index),'b.','MarkerSize',10)
 xlabel('Non-ribosomal proteins'), ylabel('J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 13;
 f3.WindowState = 'maximized';
 

f4 = figure(4); 
subplot(221)
plot(log10(a_ratio_non_ribosomal_sorted),'b.','MarkerSize',15)
 xlabel('Non-ribosomal proteins'), ylabel('Log10 J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
subplot(222)
plot(log10(a_ratio_ribosomal_sorted),'b.','MarkerSize',15)
 xlabel('Ribosomal proteins'), ylabel('Log10 J_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
%  subplot(233)
% plot(a_ratio_all_sorted,'b.','MarkerSize',10)
%  xlabel('All proteins'), ylabel('J_k(r=1) ')
%  grid on
%  ax = gca;
%  ax.FontSize = 13;
 subplot(223)
plot(log10(a_ratio_non_ribosomal_sorted./Length_proteins_non_ribosomal(a_ratio_non_ribosomal_sorted_index)),'b.','MarkerSize',15)
 xlabel('Non-ribosomal proteins'), ylabel('Log10 J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
subplot(224)
plot(log10(a_ratio_ribosomal_sorted./Length_proteins_ribosomal(a_ratio_ribosomal_sorted_index)),'b.','MarkerSize',15)
 xlabel('Ribosomal proteins'), ylabel('Log10 J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
%  subplot(236)
% plot(a_ratio_all_sorted./Length_proteins_all(a_ratio_all_sorted_index),'b.','MarkerSize',10)
%  xlabel('All proteins'), ylabel('J_k/l_k(r=1) ')
%  grid on
%  ax = gca;
%  ax.FontSize = 13;
 %print -dpng ../Texte/images/LOG_J_sorted.png
 
%  f4.WindowState = 'maximized';
%   exportgraphics(f4,'../Texte/images/LOG_J_sorted.png','Resolution',300)
  
  f5 = figure(5); 
subplot(211)
plot(log10(a_ratio_non_ribosomal_sorted./Length_proteins_non_ribosomal(a_ratio_non_ribosomal_sorted_index)),'b.','MarkerSize',15)
 xlabel('Non-ribosomal proteins'), ylabel('Log10 J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
subplot(212)
plot(log10(a_ratio_ribosomal_sorted./Length_proteins_ribosomal(a_ratio_ribosomal_sorted_index)),'b.','MarkerSize',15)
 xlabel('Ribosomal proteins'), ylabel('Log10 J_k/l_k(r=1) ')
 grid on
 ax = gca;
 ax.FontSize = 16;
%  subplot(236)
% plot(a_ratio_all_sorted./Length_proteins_all(a_ratio_all_sorted_index),'b.','MarkerSize',10)
%  xlabel('All proteins'), ylabel('J_k/l_k(r=1) ')
%  grid on
%  ax = gca;
%  ax.FontSize = 13;
 %print -dpng ../Texte/images/LOG_J_sorted.png
 
%  f5.WindowState = 'maximized';
%   exportgraphics(f5,'../Texte/images/LOG_J_sorted2.png','Resolution',300)


 %%%%%%%%%%%% AVERAGE J:
 
 %There are 3551 non-ribosomal and 68 ribosomal proteins in the data
 