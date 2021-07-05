%% Set-up
% Need to run MainAnalyses.m before this to have individual correlation coefficients

cols = cbrewer('qual', 'Set2', 3); 

%% Plot of average vs modal controllability

controlsavecontrank = [controls.connectome.avecontrank];
controlsavecontrank = controlsavecontrank(:);
controlsmodalcontrank = [controls.connectome.modalcontrank];
controlsmodalcontrank = controlsmodalcontrank(:);
controlsline = fit(controlsavecontrank,controlsmodalcontrank,'poly1');

patientsavecontrank = [patients.connectome.avecontrank];
patientsavecontrank = patientsavecontrank(:);
patientsmodalcontrank = [patients.connectome.modalcontrank];
patientsmodalcontrank = patientsmodalcontrank(:);
patientsline = fit(patientsavecontrank,patientsmodalcontrank,'poly1');

vnsavecontrank = [vns.connectome.avecontrank];
vnsavecontrank = vnsavecontrank(:);
vnsmodalcontrank = [vns.connectome.modalcontrank];
vnsmodalcontrank = vnsmodalcontrank(:);
vnsline = fit(vnsavecontrank,vnsmodalcontrank,'poly1');

subplot(3,3,[1 2])
% scatter(mean([controls.connectome.avecontrank],2),mean([controls.connectome.modalcontrank],2),[],cols(1,:),'filled')
scatter(controlsavecontrank,controlsmodalcontrank,5,cols(1,:),'filled');
hold on
line = plot(controlsline,'predfunc');
line(1).Color = cols(1,:);
line(2).Color = cols(1,:);
line(3).Color = cols(1,:);
line(1).LineWidth = 2;
line(2).LineWidth = 2;
line(3).LineWidth = 2;
legend('off');
%title('Controls')
xlim([0 253])
ylim([0 253])
xlabel('Rank Average Contollability')
ylabel('Rank Modal Controllability')
set(gca,'FontSize',15)
box off

subplot(3,3,[4 5])
% scatter(mean([patients.connectome.avecontrank],2),mean([patients.connectome.modalcontrank],2),[],cols(1,:),'filled')
scatter(patientsavecontrank,patientsmodalcontrank,5,cols(2,:),'filled');
hold on
line = plot(patientsline,'predfunc');
line(1).Color = cols(2,:);
line(2).Color = cols(2,:);
line(3).Color = cols(2,:);
line(1).LineWidth = 2;
line(2).LineWidth = 2;
line(3).LineWidth = 2;
legend('off');
%title('patients')
xlim([0 253])
ylim([0 253])
xlabel('Rank Average Contollability')
ylabel('Rank Modal Controllability')
set(gca,'FontSize',15)
box off


subplot(3,3,[7 8])
% scatter(mean([vns.connectome.avecontrank],2),mean([vns.connectome.modalcontrank],2),[],cols(1,:),'filled')
scatter(vnsavecontrank,vnsmodalcontrank,5,cols(3,:),'filled');
hold on
line = plot(vnsline,'predfunc');
line(1).Color = cols(3,:);
line(2).Color = cols(3,:);
line(3).Color = cols(3,:);
line(1).LineWidth = 2;
line(2).LineWidth = 2;
line(3).LineWidth = 2;
legend('off');
%title('vns')
xlim([0 253])
ylim([0 253])
xlabel('Rank Average Contollability')
ylabel('Rank Modal Controllability')
set(gca,'FontSize',15)
box off

subplot(3,3,3)
violin([controls.connectome.avmcorr]','mc',[],'medc',[],'facecolor',cols(1,:));
ylim([-0.9 0.3])
ylabel('A-M Correlation')
set(gca,'FontSize',15)
box off
set(gca,'Xcolor','none')
yline(0)

subplot(3,3,6)
violin([patients.connectome.avmcorr]','mc',[],'medc',[],'facecolor',cols(2,:));
ylim([-0.9 0.3])
ylabel('A-M Correlation')
set(gca,'FontSize',15)
box off
set(gca,'Xcolor','none')
yline(0)

subplot(3,3,9)
violin([vns.connectome.avmcorr]','mc',[],'medc',[],'facecolor',cols(3,:));
ylim([-0.9 0.3])
ylabel('A-M Correlation')
set(gca,'FontSize',15)
box off
set(gca,'Xcolor','none')
yline(0)

%% Save

saveas(gcf,'FinalFigs/Figure3.png')
