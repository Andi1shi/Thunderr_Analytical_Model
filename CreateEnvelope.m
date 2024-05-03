function CreateEnvelope(T,Nruns,BestFitness,BestFitIterMat,Font)
% Statistical analysis of BestFitness
Stat(1) = min(BestFitness);               % Determing the best fitness function value
Stat(2) = max(BestFitness);               % Determing the worst fitness function value
Stat(3) = mean(BestFitness);              % Determing the mean fitness function value
Stat(4) = std(BestFitness);               % Determing the standard deviation
BestFitIterMean = mean(BestFitIterMat);
BesFitIterStd = std(BestFitIterMat);
env = [min(BestFitIterMat,[],1) ; max(BestFitIterMat,[],1)];

fig = figure;
fig.WindowState = 'maximized';
iterations = 0:T;
%plot(iterations,BestFitIterMat,'.-')
hold on
p1 = plot(iterations,BestFitIterMean + BesFitIterStd,':b','LineWidth',3);
p2 = plot(iterations,BestFitIterMean,'--k','LineWidth',3);
p3 = plot(iterations,BestFitIterMean - BesFitIterStd,':g','LineWidth',3);
p4 = plot(iterations, env(1,:), '-m', 'LineWidth',3);
p5 = plot(iterations, env(2,:), '-r', 'LineWidth',3);
x2 = [iterations, fliplr(iterations)];
inBetween = [env(1,:), fliplr(env(2,:))];
fill(x2, inBetween, 'c','FaceAlpha',0.2);
grid minor
xlabel('T (Iteration Number)')
ylabel('F (Objective Function Value)')
title (['ALGORITHM: TLBO   Total Run = ' num2str(Nruns)])

str1 = ['F_{min} = ',num2str(round(Stat(1),3))];
str2 = ['F_{max} = ',num2str(round(Stat(2),3))];
str3 = ['m_F = ',num2str(round(Stat(3),3))];

str4 = ['m_F',' + ','s_F = ',num2str(round(Stat(3) + Stat(4),3))];
str5 = ['m_F',' - ','s_F = ',num2str(round(Stat(3) - Stat(4),3))];
tx = T;
tyMin = Stat(1);
tyMax = Stat(2);
tyMean = Stat(3);

tyStd = Stat(4);
ss = 15;
text(tx, tyMin - 0.025*tyMean ,str1,'Color','k','FontSize',ss,'HorizontalAlignment','left')
plot(tx,tyMin,'ok','MarkerFaceColor','k','MarkerSize',4)
text(tx, tyMax + 0.025*tyMean ,str2,'Color','k','FontSize',ss,'HorizontalAlignment','left')
plot(tx,tyMax,'ok','MarkerFaceColor','k','MarkerSize',4)
text(tx, tyMean,str3,'Color','k','FontSize',ss,'HorizontalAlignment','left')
plot(tx,tyMean,'ok','MarkerFaceColor','k','MarkerSize',4)
text(tx, tyMean + tyStd - 0.025*tyMean ,str4,'Color','k','FontSize',ss,'HorizontalAlignment','left')
plot(tx,tyMean + tyStd,'ok','MarkerFaceColor','k','MarkerSize',4)
text(tx, tyMean - tyStd + 0.025*tyMean,str5,'Color','k','FontSize',ss,'HorizontalAlignment','left')
plot(tx,tyMean - tyStd,'ok','MarkerFaceColor','k','MarkerSize',4)
legend([p5,p1,p2,p3,p4],{'Envelope Upper Bound','Mean + Std','Mean','Mean - Std','Envelope Lower Bound'},'FontSize',Font)
ax = gca;
ax.FontSize = Font;
