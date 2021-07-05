%% Fluff and setup
% to be run after MainAnalysis.m

cd /Controllability;  % input location of data structures
cols = cbrewer('qual', 'Paired', 6);     % color scheme
Np  = size(patients.connectome,2);                % Number of patients
jit = (rand(Np,1)-0.5)/4;              % jitter for dot plots
outcomes = [patients.connectome.Engel] == 1; 

postop = load('PostopPatients/controllabilities.mat');

%% Figure5: To visualise resections. 
% Make resection structure, assume control resections are 0?

for a = 1:length(patients.connectome)
    resection.zwdeg(a) = patients.connectome(a).resection.zwdeg
    resection.zavecont(a) = patients.connectome(a).resection.zavecont
    resection.zmodalcont(a) = patients.connectome(a).resection.zmodalcont
    
    resection.zwdegrob(a) = patients.connectome(a).resection.zwdegrob
    resection.zavecontrob(a) = patients.connectome(a).resection.zavecontrob
    resection.zmodalcontrob(a) = patients.connectome(a).resection.zmodalcontrob
    
    resection.zwdegwb(a) = patients.connectome(a).resection.zwdegwb
    resection.zavecontwb(a) = patients.connectome(a).resection.zavecontwb 
    resection.zmodalcontwb(a) = patients.connectome(a).resection.zmodalcontwb
end

%% Z-score post-op parcels

for a = 1:length(postop.connectome)
    
    % wdegrank
    
    [~,r] = sort(postop.connectome(a).wdeg);
    rankV(r) = 1:numel(postop.connectome(a).wdeg);
    postop.connectome(a).wdegrank = (rankV)';
    
    % avecontrank
    
    [~,r] = sort(postop.connectome(a).avecont);
    rankV(r) = 1:numel(postop.connectome(a).avecont);
    postop.connectome(a).avecontrank = (rankV)';
    
    % modalcontrank
   
    [~,r] = sort(postop.connectome(a).modalcont);
    rankV(r) = 1:numel(postop.connectome(a).modalcont);
    postop.connectome(a).modalcontrank = (rankV)';
    
end

for a = 1:length(postop.connectome)
    for b = 1:253
        zwdeg(b) = (postop.connectome(a).wdegrank(b) - controlwdegrankmean(b))/controlwdegrankstd(b);
        zavecont(b) = (postop.connectome(a).avecontrank(b) - controlavecontrankmean(b))/controlavecontrankstd(b);
        zmodalcont(b) = (postop.connectome(a).modalcontrank(b) - controlmodalcontrankmean(b))/controlmodalcontrankstd(b); 
    end
    
    postop.connectome(a).zwdeg = zwdeg';
    postop.connectome(a).zavecont = zavecont';
    postop.connectome(a).zmodalcont = zmodalcont';
    
    clear zwdeg zavecont zmodalcont;
    
end

%% Calculate Metrics: For post-op, it will be only for remaining parcels

for a = 1:length(postop.connectome)
    
    resection.zwdegpo(a) = mean(postop.connectome(a).zwdeg);
    resection.zavecontpo(a) = mean(postop.connectome(a).zavecont);
    resection.zmodalcontpo(a) = mean(postop.connectome(a).zmodalcont);
    
end

%% Effect of resected parcels

subplot(5,3,1)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zwdegwb(a) resection.zwdeg(a)],'color',cols(1+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zwdegwb, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zwdeg, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')

errorbar([0.6],mean(resection.zwdegwb(outcomes==0)),1.96*std(resection.zwdegwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zwdegwb(outcomes==1)),1.96*std(resection.zwdegwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zwdeg(outcomes==0)),1.96*std(resection.zwdeg(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zwdeg(outcomes==1)),1.96*std(resection.zwdeg(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    

xlim([.3 2.7])
ylim([-3.1 3.1])
yticks([-2 -1 0 1 2])
yticklabels({'-2.0','-1.0','0','1.0','2.0'})
sigstar([.8 2.4])
yline(0)
set(gca,'xticklabel',[])
xticks([])
ylabel('Z-Score')
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

subplot(5,3,2)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zavecontwb(a) resection.zavecont(a)],'color',cols(3+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zavecontwb, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zavecont, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
errorbar([0.6],mean(resection.zavecontwb(outcomes==0)),1.96*std(resection.zavecontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zavecontwb(outcomes==1)),1.96*std(resection.zavecontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zavecont(outcomes==0)),1.96*std(resection.zavecont(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zavecont(outcomes==1)),1.96*std(resection.zavecont(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    
xlim([.3 2.7])
ylim([-3.1 3.1])
yticks([-2 -1 0 1 2])
yticklabels({'-2.0','-1.0','0','1.0','2.0'})
yline(0)
set(gca,'xticklabel',[])
xticks([])
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

subplot(5,3,3)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zmodalcontwb(a) resection.zmodalcont(a)],'color',cols(5+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zmodalcontwb, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zmodalcont, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
errorbar([0.6],mean(resection.zmodalcontwb(outcomes==0)),1.96*std(resection.zmodalcontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zmodalcontwb(outcomes==1)),1.96*std(resection.zmodalcontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zmodalcont(outcomes==0)),1.96*std(resection.zmodalcont(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zmodalcont(outcomes==1)),1.96*std(resection.zmodalcont(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    
xlim([.3 2.7])
ylim([-3.1 3.1])
yticks([-2 -1 0 1 2])
yticklabels({'-2.0','-1.0','0','1.0','2.0'})
sigstar([.8 2.4])
yline(0)
set(gca,'xticklabel',[])
xticks([])
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

%% Effect of on rest of brain

subplot(5,3,4)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zwdegwb(a) resection.zwdegrob(a)],'color',cols(1+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zwdegwb, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zwdegrob, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')
errorbar([0.6],mean(resection.zwdegwb(outcomes==0)),1.96*std(resection.zwdegwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zwdegwb(outcomes==1)),1.96*std(resection.zwdegwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zwdegrob(outcomes==0)),1.96*std(resection.zwdegrob(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zwdegrob(outcomes==1)),1.96*std(resection.zwdegrob(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    
xlim([.3 2.7])
ylim([-.5 .5])
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
yline(0)
sigstar([.8 2.4])
set(gca,'xticklabel',[])
xticks([])
ylabel('Z-Score')
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

subplot(5,3,5)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zavecontwb(a) resection.zavecontrob(a)],'color',cols(3+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zavecontwb, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zavecontrob, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
errorbar([0.6],mean(resection.zavecontwb(outcomes==0)),1.96*std(resection.zavecontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zavecontwb(outcomes==1)),1.96*std(resection.zavecontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zavecontrob(outcomes==0)),1.96*std(resection.zavecontrob(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zavecontrob(outcomes==1)),1.96*std(resection.zavecontrob(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    
xlim([.3 2.7])
ylim([-.5 .5])
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
yline(0)
set(gca,'xticklabel',[])
xticks([])
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

subplot(5,3,6)
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zmodalcontwb(a) resection.zmodalcontrob(a)],'color',cols(5+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zmodalcontwb, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zmodalcontrob, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
errorbar([0.6],mean(resection.zmodalcontwb(outcomes==0)),1.96*std(resection.zmodalcontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zmodalcontwb(outcomes==1)),1.96*std(resection.zmodalcontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zmodalcontrob(outcomes==0)),1.96*std(resection.zmodalcontrob(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zmodalcontrob(outcomes==1)),1.96*std(resection.zmodalcontrob(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    
xlim([.3 2.7])
ylim([-.5 .5])
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
yline(0)
sigstar([.8 2.4])
set(gca,'xticklabel',[])
xticks([])
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

%% Postop connectime with TCKEXCLUDE

subplot(5,3,7)      % Weighted Degree
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.zwdegwb(a) resection.zwdegpo(a)],'color',cols(1+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.zwdegwb, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.zwdegpo, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')

errorbar([0.6],mean(resection.zwdegwb(outcomes==0)),1.96*std(resection.zwdegwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([0.8],mean(resection.zwdegwb(outcomes==1)),1.96*std(resection.zwdegwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    
errorbar([2.2],mean(resection.zwdegpo(outcomes==0)),1.96*std(resection.zwdegpo(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(1,:),'LineWidth',4)    
errorbar([2.4],mean(resection.zwdegpo(outcomes==1)),1.96*std(resection.zwdegpo(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(2,:),'LineWidth',4)    

xlim([.3 2.7])
ylim([-.5 .5])
sigstar({[.8 2.4],[.6 2.2]})
yline(0)
set(gca,'xticklabel',[])
xticks([])
ylabel('Z-Score')
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')

subplot(5,3,8)      % Average Controllability
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.avecontwb(a) resection.avecontpo(a)],'color',cols(3+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.avecontwb, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.avecontpo, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')

errorbar([0.6],mean(resection.avecontwb(outcomes==0)),1.96*std(resection.avecontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([0.8],mean(resection.avecontwb(outcomes==1)),1.96*std(resection.avecontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    
errorbar([2.2],mean(resection.avecontpo(outcomes==0)),1.96*std(resection.avecontpo(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(3,:),'LineWidth',4)    
errorbar([2.4],mean(resection.avecontpo(outcomes==1)),1.96*std(resection.avecontpo(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(4,:),'LineWidth',4)    

xlim([.3 2.7])
ylim([-.5 .5])
sigstar({[.8 2.4]})
yline(0)
set(gca,'xticklabel',[])
xticks([])
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')


subplot(5,3,9)         % Modal Controllability
hold on
for a=1:52
    plot([1.2+jit(a) 1.8+jit(a)],[resection.modalcontwb(a) resection.modalcontpo(a)],'color',cols(5+outcomes(a),:))
end
scatter((ones(Np,1).*1.2)+jit, resection.modalcontwb, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
scatter((ones(Np,1).*1.8)+jit, resection.modalcontpo, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')

errorbar([0.6],mean(resection.modalcontwb(outcomes==0)),1.96*std(resection.modalcontwb(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([0.8],mean(resection.modalcontwb(outcomes==1)),1.96*std(resection.modalcontwb(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    
errorbar([2.2],mean(resection.modalcontpo(outcomes==0)),1.96*std(resection.modalcontpo(outcomes==0))/sqrt(19),'+','MarkerSize',10,'Color',cols(5,:),'LineWidth',4)    
errorbar([2.4],mean(resection.modalcontpo(outcomes==1)),1.96*std(resection.modalcontpo(outcomes==1))/sqrt(33),'+','MarkerSize',10,'Color',cols(6,:),'LineWidth',4)    

xlim([.3 2.7])
ylim([-.5 .5])
sigstar({[.8 2.4],[.6 2.2]})
yline(0)
set(gca,'xticklabel',[])
xticks([])
yticks([-.5 -.25 0 .25 .5])
yticklabels({'-.50','-.25','0','.25','.50'})
set(gca,'FontSize',13)
%set(gca,'Ydir','reverse')
set(gca,'XColor', 'none')


%% Virtual resection

for a = 1:length(patients.connectome)
    resection.zofzwdeg(a) = patients.connectome(a).resection.zofzwdeg;
    resection.zofzavecont(a) = patients.connectome(a).resection.zofzavecont;
    resection.zofzmodalcont(a) = patients.connectome(a).resection.zofzmodalcont;
end

    subplot(5,3,[10 13])                      % WDEG
    hold on
    % plot([0.8, 1.2], [mean(resection.zofzwdeg(outcomes==1)), mean(resection.zofzwdeg(outcomes==1))], 'linewidth', 4, 'color', cols(2,:))
    % plot([0.8, 1.2], [mean(resection.zofzwdeg(outcomes==0)), mean(resection.zofzwdeg(outcomes==0))], 'linewidth', 4, 'color', cols(1,:))
    errorbar([1.3],mean(resection.zofzwdeg(outcomes==1)),1.96*std(resection.zofzwdeg(outcomes==1))/sqrt(33),'+','MarkerSize',20,'Color',cols(2,:),'LineWidth',4)    
    errorbar([1.5],mean(resection.zofzwdeg(outcomes==0)),1.96*std(resection.zofzwdeg(outcomes==0))/sqrt(19),'+','MarkerSize',20,'Color',cols(1,:),'LineWidth',4)    
    scatter(ones(Np,1)+jit, resection.zofzwdeg, [], cols((1)+outcomes,:), 'filled','MarkerEdgeColor','k')
    ylim([-7 7])
    xlim([.8 1.6])
    yline(0)
  
    %legend('Engel Class I','Engel Class II-IV','location','southeast')
    ylabel('Z-Score compared to 1000 virtual resections')
    set(gca,'xticklabel',[])
    set(gca,'XColor', 'none')
    %set(gca,'Ydir','reverse')
    set(gca,'FontSize',13)
    
    subplot(5,3,[11 14])                      % avecont
    hold on
    % plot([0.8, 1.2], [mean(resection.zofzavecont(outcomes==1)), mean(resection.zofzavecont(outcomes==1))], 'linewidth', 4, 'color', cols(2,:))
    % plot([0.8, 1.2], [mean(resection.zofzavecont(outcomes==0)), mean(resection.zofzavecont(outcomes==0))], 'linewidth', 4, 'color', cols(1,:))
    errorbar([1.3],mean(resection.zofzavecont(outcomes==1)),1.96*std(resection.zofzavecont(outcomes==1))/sqrt(33),'+','MarkerSize',20,'Color',cols(4,:),'LineWidth',4)    
    errorbar([1.5],mean(resection.zofzavecont(outcomes==0)),1.96*std(resection.zofzavecont(outcomes==0))/sqrt(19),'+','MarkerSize',20,'Color',cols(3,:),'LineWidth',4)    
    scatter(ones(Np,1)+jit, resection.zofzavecont, [], cols((3)+outcomes,:), 'filled','MarkerEdgeColor','k')
    ylim([-7 7])
    xlim([.8 1.6])
    yline(0)

    %legend('Engel Class I','Engel Class II-IV','location','southeast')
    set(gca,'xticklabel',[])
    set(gca,'XColor', 'none')
   % set(gca,'Ydir','reverse')
    set(gca,'FontSize',13)
    
    subplot(5,3,[12 15])                      % modalcont
    hold on
    % plot([0.8, 1.2], [mean(resection.zofzmodalcont(outcomes==1)), mean(resection.zofzmodalcont(outcomes==1))], 'linewidth', 4, 'color', cols(2,:))
    % plot([0.8, 1.2], [mean(resection.zofzmodalcont(outcomes==0)), mean(resection.zofzmodalcont(outcomes==0))], 'linewidth', 4, 'color', cols(1,:))
    errorbar([1.3],mean(resection.zofzmodalcont(outcomes==1)),1.96*std(resection.zofzmodalcont(outcomes==1))/sqrt(33),'+','MarkerSize',20,'Color',cols(6,:),'LineWidth',4)    
    errorbar([1.5],mean(resection.zofzmodalcont(outcomes==0)),1.96*std(resection.zofzmodalcont(outcomes==0))/sqrt(19),'+','MarkerSize',20,'Color',cols(5,:),'LineWidth',4)    
    scatter(ones(Np,1)+jit, resection.zofzmodalcont, [], cols((5)+outcomes,:), 'filled','MarkerEdgeColor','k')
    ylim([-7 7])
    xlim([.8 1.6])
    yline(0)
   
    %legend('Engel Class I','Engel Class II-IV','location','southeast')
    set(gca,'xticklabel',[])
    set(gca,'XColor', 'none')
  %  set(gca,'Ydir','reverse')
    set(gca,'FontSize',13)

%% Save

saveas(gcf,'FinalFigs/Figure6.png')
