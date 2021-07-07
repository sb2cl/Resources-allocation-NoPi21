
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               
%
% Next we obtain the relationship between free ribosomes and the 
% translation efficiency per mRNA.
% We consider that the ribosomes density decreases with the protein
% length
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

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
 large_dm= log(2)/7.5*60;
 for i=1:x
    if isnan(Dades(i,5))  
       Dades(i,5)=log10(mean_dm);
    end
 end
% An alternative considering that no mRNA half-life is larger than 7.5 mins 
% for i=1:x
%    if isnan(Dades(i,5))
%        Dades(i,5)=log10(general_dm);
%    elseif (Dades(i,5)<log10(large_dm))
%        Dades(i,5)= log10(large_dm);
%    end
% end

% We consider ribosomal proteins independently
Ribosomal_indices = 2058:2125;

% Get rid of genes at which some of the values bp,bm,lp is NaN
Bad_genes=any(isnan(Dades(:,1:3)),2);
Good_genes_all=find(1-Bad_genes); %Gives their indices 
tempo=Dades(Good_genes_all,:);

Good_genes_ribosomal=Good_genes_all(ismember(Good_genes_all,Ribosomal_indices));
Good_genes_non_ribosomal = Good_genes_all(not(ismember(Good_genes_all,Ribosomal_indices')));

% Transcription, translation rates and protein length for all goodies:
Translation_rates_all=10.^Dades(Good_genes_all,1)/60;    %protein/(mRNA*min)
Length_proteins_all= 10.^Dades(Good_genes_all,3);         % amino acids (aa)
Transcription_rates_all=10.^Dades(Good_genes_all,2)/60;  %mRNA/min
Transcription_rates_all_aa= Transcription_rates_all.*Length_proteins_all; %aa/min
R_density_all=1./(14.223*Length_proteins_all.^0.097);
Mean_R_density_all=mean(R_density_all)

% Transcription, translation rates and protein length for all ribosomal goodies:
Translation_rates_ribosomal=10.^Dades(Good_genes_ribosomal,1)/60;    %protein/(mRNA*min)
Length_proteins_ribosomal= 10.^Dades(Good_genes_ribosomal,3);         % amino acids (aa)
Transcription_rates_ribosomal=10.^Dades(Good_genes_ribosomal,2)/60;  %mRNA/min
Transcription_rates_ribosomal_aa= Transcription_rates_ribosomal.*Length_proteins_ribosomal; %aa/min
Names_ribosomes= Names_genes(Good_genes_ribosomal);
Total_length_proteins_ribosomal=sum(Length_proteins_ribosomal);
%R_density_ribosomal=4.5e-2-0.1e-4*(Length_proteins_ribosomal-100);
R_density_ribosomal=1./(14.223*Length_proteins_ribosomal.^0.097);
Mean_R_density_ribosomal=mean(R_density_ribosomal);

% Transcription, translation rates and protein length for all NON ribosomal goodies:
Translation_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,1)/60;    %protein/(mRNA*min)
Length_proteins_non_ribosomal= 10.^Dades(Good_genes_non_ribosomal,3);         % amino acids (aa)
Transcription_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,2)/60;  %mRNA/min
Transcription_rates_non_ribosomal_aa= Transcription_rates_non_ribosomal.*Length_proteins_non_ribosomal; %aa/min
%R_density_non_ribosomal=4.5e-2-0.1e-4*(Length_proteins_non_ribosomal-100);
R_density_non_ribosomal=1./(14.223*Length_proteins_non_ribosomal.^0.097);
Mean_R_density_non_ribosomal=mean(R_density_non_ribosomal);

% Mean mRNA degradation rate
mRNA_degradation_rates_all=10.^Dades(Good_genes_all,5)/60;  %1/min
mRNA_degradation_rates_ribosomal=10.^Dades(Good_genes_ribosomal,5)/60;  %1/min
mRNA_degradation_rates_non_ribosomal=10.^Dades(Good_genes_non_ribosomal,5)/60;  %1/min

%%%%%%% GENERAL CELL DATA:

mu_Hausser= log(2)/21.5; %0.032
maa = 182.6e-24; %average weight for an amino acid
% For the whole cell, Hausser 2019 estimates:
mc_Hausser = 180e-15;
nu=1260;
b=maa*nu/mc_Hausser; %b=5.3506e-07

td=20:5:100;
mu = log(2)./td;
Phib_ra = mu/b;
% Range of r assuming N(350,80) -4sigma to +4sigma
% R.J. Harvey, ﻿Fraction of ribosomes synthesizing protein as a function of
% specific growth rate. J. ﻿of Bacteriology, 1973. 
free_r = 30:1:670;
RBS_strength= 0.005:0.005:0.3;

%%%%%%% FOR HAUSSER DATA:
[X,Y] = meshgrid(free_r,RBS_strength);

YpmRNA_Hausser_non_ribosomal = Translation_rates_non_ribosomal./mRNA_degradation_rates_non_ribosomal;
Z_non_ribosomal=[];
for k=1:size(Translation_rates_non_ribosomal,1)
    le= 1/R_density_non_ribosomal(k);  
    Z_non_ribosomal(k).z_values = 0.62*nu/le*X./(mRNA_degradation_rates_non_ribosomal(k)./Y + mu_Hausser*X);
end
Mean_YpmRNA_Hausser_non_ribosomal=mean(YpmRNA_Hausser_non_ribosomal);
Mean_Z_non_ribosomal.z_values =  0.62*nu*Mean_R_density_non_ribosomal*X./(mean(mRNA_degradation_rates_non_ribosomal)./Y + mu_Hausser*X);
Mean_Z_non_ribosomal.z_values_approx =  0.62*nu*Mean_R_density_non_ribosomal*X./(mean(mRNA_degradation_rates_non_ribosomal)./Y + 0*mu_Hausser*X);

YpmRNA_Hausser_ribosomal = Translation_rates_ribosomal./mRNA_degradation_rates_ribosomal;
Z_ribosomal=[];
for k=1:size(Translation_rates_ribosomal,1)
    le= 1/R_density_ribosomal(k); 
    Z_ribosomal(k).z_values = 0.62*nu/le*X./(mRNA_degradation_rates_ribosomal(k)./Y + mu_Hausser*X);
end
Mean_YpmRNA_Hausser_ribosomal=mean(YpmRNA_Hausser_ribosomal);
Mean_Z_ribosomal.z_values = 0.62*nu*Mean_R_density_ribosomal*X./(mean(mRNA_degradation_rates_ribosomal)./Y + mu_Hausser*X);
Mean_Z_ribosomal.z_values_approx = 0.62*nu*Mean_R_density_ribosomal*X./(mean(mRNA_degradation_rates_ribosomal)./Y + 0*mu_Hausser*X);

YpmRNA_Hausser_all = Translation_rates_all./mRNA_degradation_rates_all;
Z_all=[];
for k=1:size(Translation_rates_all,1)
    le= 1/R_density_all(k);  
    Z_all(k).z_values = 0.62*nu/le*X./(mRNA_degradation_rates_all(k)./Y + mu_Hausser*X);
end
Mean_YpmRNA_Hausser_all=mean(YpmRNA_Hausser_all);
Mean_Z_all.z_values = 0.62*nu*Mean_R_density_all*X./(mean(mRNA_degradation_rates_all)./Y + 0*mu_Hausser*X);

ColorMap = [
    153,0,0;
    255,0,0;
    255,102,102;
    255, 204,204]./255;

Fontsize=16;

f1=figure(1)
subplot(121)
M_non_ribosomal=[];
c_non_ribosoma=[];
for k=1:size(Translation_rates_non_ribosomal,1)
    % plot for r_density=3.5 ribosomes per 100 codons (Bionumbers)
    [M_non_ribosomal(k).data,c_non_ribosomal(k).data]= contour(X,Y,Z_non_ribosomal(k).z_values,[YpmRNA_Hausser_non_ribosomal(k),YpmRNA_Hausser_non_ribosomal(k)]);
    hold on
end
 [M_mean_non_ribosomal,c_mean_non_ribosomal]= contour(X,Y,Mean_Z_non_ribosomal.z_values,[Mean_YpmRNA_Hausser_non_ribosomal,Mean_YpmRNA_Hausser_non_ribosomal],'Linewidth',3,'Color','r');
 [M_mean_non_ribosomal_approx,c_mean_non_ribosomal_approx]=contour(X,Y,Mean_Z_non_ribosomal.z_values_approx,[Mean_YpmRNA_Hausser_non_ribosomal,Mean_YpmRNA_Hausser_non_ribosomal],'Linewidth',3,'Color','g');

 xlabel('r, Non-ribosomal','Interpreter','latex','FontSize',Fontsize), ylabel('$K_C^0$','Interpreter','latex','FontSize',Fontsize)
 %title('K_{p/mRNA} (molec. min^{-1} mRNA^{-1})')
 grid on
 ax = gca;
 ax.FontSize = 16;
 ax.YLim = [0.0 0.25];
 hold off
 
 subplot(122)
 M_ribosomal=[];
 c_ribosoma=[];
for k=1:size(Translation_rates_ribosomal,1)
    % plot for le=30 codons
    [M_ribosomal(k).data,c_ribosomal(k).data]= contour(X,Y,Z_ribosomal(k).z_values,[YpmRNA_Hausser_ribosomal(k),YpmRNA_Hausser_ribosomal(k)]);
    hold on
end
 [M_mean_ribosomal,c_mean_ribosomal]= contour(X,Y,Mean_Z_ribosomal.z_values,[Mean_YpmRNA_Hausser_ribosomal,Mean_YpmRNA_Hausser_ribosomal],'Linewidth',3,'Color','r');
 [M_mean_ribosomal_approx,c_mean_ribosomal_approx]= contour(X,Y,Mean_Z_ribosomal.z_values_approx,[Mean_YpmRNA_Hausser_ribosomal,Mean_YpmRNA_Hausser_ribosomal],'Linewidth',3,'Color','g');
 xlabel('r, Ribosomal','Interpreter','latex','FontSize',Fontsize), ylabel('$K_C^0$','Interpreter','latex','FontSize',Fontsize)
 %title('K_{p/mRNA} (molec. min^{-1} mRNA^{-1})')
 grid on
 ax = gca;
 ax.FontSize = 16;
 ax.YLim = [0.0 0.25];
 hold off
 %print -dpng ../Texte/images/f_effective_translation_r_Ypmrna_density2.png
% f1.WindowState = 'maximized';
%exportgraphics(f1,'../images_def/f_effective_translation_r_Ypmrna_density2.png','Resolution',300)
 
 %%%%% EXTRACT INDIVIDUAL CURVES AND LOG PLOT FOR X-AXIS (FREE RIBOSOMES):
 
 
 %%% LOG PLOTS %%%%%
 
 %M_non_ribosomal
 
 figure(2)
 subplot(121)
 for k=1:length(M_non_ribosomal)
   plot(log10(M_non_ribosomal(k).data(1,2:end)),log10(M_non_ribosomal(k).data(2,2:end)),'Linewidth',1)
   hold on
 end
 xlabel('$\log_{10}r$','Interpreter','latex','FontSize',Fontsize), ylabel('$\log_{10}K_C^0$, non-ribosomal','Interpreter','latex','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 16;
 ax.YLim = [-2 -0.5];
 hold off
 
 subplot(122)
 for k=1:length(M_ribosomal)
   plot(log10(M_ribosomal(k).data(1,2:end)),log10(M_ribosomal(k).data(2,2:end)),'Linewidth',1)
   hold on
 end
  xlabel('$\log_{10}r$','Interpreter','latex','FontSize',Fontsize),  ylabel('$\log_{10}K_C^0$, ribosomal','Interpreter','latex','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 16;
 ax.YLim = [-2 -0.5];
 hold off
 
 figure(3)
 subplot(121)
 for k=1:length(M_non_ribosomal)
   plot(log10(M_non_ribosomal(k).data(1,2:end)),M_non_ribosomal(k).data(2,2:end),'Linewidth',1)
   hold on
 end
 xlabel('$r$','Interpreter','latex','FontSize',Fontsize), ylabel('$\log_{10}K_C^0$, non-ribosomal','Interpreter','latex','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 16;
% ax.YLim = [-2 -0.5];
 hold off
 
 subplot(122)
 for k=1:length(M_ribosomal)
   plot(log10(M_ribosomal(k).data(1,2:end)),M_ribosomal(k).data(2,2:end),'Linewidth',1)
   hold on
 end
 xlabel('$r$','Interpreter','latex','FontSize',Fontsize), ylabel('$\log_{10}K_C^0$, ribosomal','Interpreter','latex','FontSize',Fontsize)
 grid on
 ax = gca;
 ax.FontSize = 16;
 %ax.YLim = [-2 -0.5];
 hold off
 
 
