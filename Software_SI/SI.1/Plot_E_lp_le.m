

k=0.62;

f1=figure(1)
le= 15;
hand1=fplot(@(x)x/le*( 1 - (x./(x+le)).^(x/le) ),[200 500],'Linewidth',3,'Color','b');
hold on, grid on
hand2=fplot(@(x)x/le*k,[200 500],'-.','Linewidth',3,'Color','k');

le= 20;
hand3=fplot(@(x)x/le*( 1 - (x./(x+le)).^(x/le) ),[200 500],'Linewidth',3,'Color','g');
hand4=fplot(@(x)x/le*k,[200 500],'-.','Linewidth',3,'Color','k');

le= 25;
hand5=fplot(@(x)x/le*( 1 - (x./(x+le)).^(x/le) ),[200 500],'Linewidth',3,'Color','y');
hand6=fplot(@(x)x/le*k,[200 500],'-.','Linewidth',3,'Color','k');

le= 30;
hand7=fplot(@(x)x/le*( 1 - (x./(x+le)).^(x/le) ),[200 500],'Linewidth',3,'Color','r');
hand8=fplot(@(x)x/le*k,[200 500],'-.','Linewidth',3,'Color','k');

subset=[hand1,hand3,hand5,hand7]; % This way only legends for the subset of plots in handles.

legend(subset,'$l_e=15$','$l_e=20$','$l_e=25$','$l_e=30$','Location','northwest','FontSize',16,'Interpreter','latex','AutoUpdate','off')
xlabel('$l_{pk}$','FontSize',18,'Interpreter','latex'), ylabel('$E_{mk}(l_{pk},l_e)$','FontSize',18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
hold off

%exportgraphics(f1,'../images_def/Emax_lxle.png','Resolution',300)

 
