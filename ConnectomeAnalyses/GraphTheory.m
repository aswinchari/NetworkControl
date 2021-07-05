%% Graph Theory Metrics: Calculates graph theory metrics to link them to the correlation coefficients
% Need to run MainAnalyses.m before this to have individual coefficients
%% Load Connectomes

patients = load('Patients/controllabilities.mat');
controls = load('Controls/controllabilities.mat'); 
vns = load('VNS/controllabilities.mat');
load('labels.mat');

%% Calculate global metrics for each graph

for a = 1:length(controls.connectome)
    
   [module,modularity] = community_louvain(controls.connectome(a).connectome);
   gefficiency = efficiency_wei(controls.connectome(a).connectome);
   [defficiency,~] = diffusion_efficiency(controls.connectome(a).connectome);
   partcoef = participation_coef(controls.connectome(a).connectome,module,0);
   density = density_und(controls.connectome(a).connectome);
   
   controls.connectome(a).module = module;
   controls.connectome(a).modularity = modularity;
   controls.connectome(a).gefficiency = gefficiency;
   controls.connectome(a).defficiency = defficiency;
   controls.connectome(a).partcoef = partcoef;
   controls.connectome(a).density = density;
   binary.controls{a} = controls.connectome(a).connectome~=0;
   
   clear modularity gefficiency defficiency module partcoef density
   
end

for a = 1:length(patients.connectome)
    
   [module,modularity] = community_louvain(patients.connectome(a).connectome);
   gefficiency = efficiency_wei(patients.connectome(a).connectome);
   [defficiency,~] = diffusion_efficiency(patients.connectome(a).connectome);
   partcoef = participation_coef(patients.connectome(a).connectome,module,0);
   density = density_und(patients.connectome(a).connectome);
   
   patients.connectome(a).module = module;
   patients.connectome(a).modularity = modularity;
   patients.connectome(a).gefficiency = gefficiency;
   patients.connectome(a).defficiency = defficiency;
   patients.connectome(a).partcoef = partcoef;
   patients.connectome(a).density = density;
   binary.patients{a} = patients.connectome(a).connectome~=0;

   
   clear modularity gefficiency defficiency module partcoef density
   
end

for a = 1:length(vns.connectome)
    
   [module,modularity] = community_louvain(vns.connectome(a).connectome);
   gefficiency = efficiency_wei(vns.connectome(a).connectome);
   [defficiency,~] = diffusion_efficiency(vns.connectome(a).connectome);
   partcoef = participation_coef(vns.connectome(a).connectome,module,0);
   density = density_und(vns.connectome(a).connectome);
   
   vns.connectome(a).module = module;
   vns.connectome(a).modularity = modularity;
   vns.connectome(a).gefficiency = gefficiency;
   vns.connectome(a).defficiency = defficiency;
   vns.connectome(a).partcoef = partcoef;
   vns.connectome(a).density = density;
   binary.vns{a} = vns.connectome(a).connectome~=0;

   
   clear modularity gefficiency defficiency module partcoef density
   
end

%% Put them all into a single structure

% Cont metrics
metrics.meanwdeg = [[controls.connectome.meanwdeg] [patients.connectome.meanwdeg] [vns.connectome.meanwdeg]]';
metrics.meanavecont = [[controls.connectome.meanavecont] [patients.connectome.meanavecont] [vns.connectome.meanavecont]]';
metrics.meanmodalcont = [[controls.connectome.meanmodalcont] [patients.connectome.meanmodalcont] [vns.connectome.meanmodalcont]]';
metrics.density = [[controls.connectome.density] [patients.connectome.density] [vns.connectome.density]]';
metrics.avecontcorr = [[controls.connectome.avecontcorr] [patients.connectome.avecontcorr] [vns.connectome.avecontcorr]]';
metrics.modalcontcorr = [[controls.connectome.modalcontcorr] [patients.connectome.modalcontcorr] [vns.connectome.modalcontcorr]]';
metrics.avmcorr = [[controls.connectome.avmcorr] [patients.connectome.avmcorr] [vns.connectome.avmcorr]]';
metrics.modularity = [[controls.connectome.modularity] [patients.connectome.modularity] [vns.connectome.modularity]]';
metrics.gefficiency = [[controls.connectome.gefficiency] [patients.connectome.gefficiency] [vns.connectome.gefficiency]]';
metrics.defficiency = [[controls.connectome.defficiency] [patients.connectome.defficiency] [vns.connectome.defficiency]]';
metrics.group = vertcat(ones(16,1), ones(52,1).*2, ones(27,1).*3);
metrics.module = [[controls.connectome.module] [patients.connectome.module] [vns.connectome.module]];
metrics.cog = [[controls.connectome.cog] [patients.connectome.cog] [vns.connectome.cog]]';
metrics.age = [[controls.connectome.age] [patients.connectome.age] [vns.connectome.age]]';


for a=1:length(patients.connectome)
patients.connectome(a).outcomes = patients.connectome(a).Engel==1
end

for a=1:length(vns.connectome)
vns.connectome(a).outcomes = vns.connectome(a).outcome==1
end

metrics.outcomes = vertcat(ones(16,1), ([patients.connectome.outcomes]+2)', ([vns.connectome.outcomes]+4)')

%% Coeff of metrics in figure below

[GTCoeff.Mod.R, GTCoeff.Mod.P, GTCoeff.Mod.RL, GTCoeff.Mod.RU] = corrcoef(metrics.avmcorr,metrics.modularity);
[GTCoeff.Geff.R, GTCoeff.Geff.P, GTCoeff.Geff.RL, GTCoeff.Geff.RU] = corrcoef(metrics.avmcorr,metrics.gefficiency);
[GTCoeff.Deff.R, GTCoeff.Deff.P, GTCoeff.Deff.RL, GTCoeff.Deff.RU] = corrcoef(metrics.avmcorr,metrics.defficiency);



%% Plot some of these against each other (Supplemental Figure)

cols = cbrewer('qual', 'Set2', 3); 
colormap(cols);

% AC vs GT metrics

subplot(4,5,1)
scatter(metrics.meanavecont,metrics.meanwdeg,50,metrics.group,'filled')
hold on
a = fit(metrics.meanavecont,metrics.meanwdeg,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
title('Mean Average Controllability')
ylabel('Mean Weighted Degree')
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,6)
scatter(metrics.meanavecont,metrics.modularity,50,metrics.group,'filled')
hold on
a = fit(metrics.meanavecont,metrics.modularity,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel([])
ylabel('Modularity')

subplot(4,5,11)
scatter(metrics.meanavecont,metrics.gefficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.meanavecont,metrics.gefficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel([])
ylabel('Global Efficiency')

subplot(4,5,16)
scatter(metrics.meanavecont,metrics.defficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.meanavecont,metrics.defficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel('Mean Average Controllability')
ylabel('Diffusion Efficiency')

% MC vs GT metrics

subplot(4,5,2)
scatter(metrics.meanmodalcont,metrics.meanwdeg,50,metrics.group,'filled')
title('Mean Modal Controllability')
hold on
a = fit(metrics.meanmodalcont,metrics.meanwdeg,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,7)
scatter(metrics.meanmodalcont,metrics.modularity,50,metrics.group,'filled')
hold on
a = fit(metrics.meanmodalcont,metrics.modularity,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,12)
scatter(metrics.meanmodalcont,metrics.gefficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.meanmodalcont,metrics.gefficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,17)
scatter(metrics.meanmodalcont,metrics.defficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.meanmodalcont,metrics.defficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel('Mean Modal Controllability')

% ACR vs GT metrics

subplot(4,5,3)
scatter(metrics.avecontcorr,metrics.meanwdeg,50,metrics.group,'filled')
hold on
a = fit(metrics.avecontcorr,metrics.meanwdeg,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
title('AC v WD Correlation')
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,8)
scatter(metrics.avecontcorr,metrics.modularity,50,metrics.group,'filled')
hold on
a = fit(metrics.avecontcorr,metrics.modularity,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,13)
scatter(metrics.avecontcorr,metrics.gefficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.avecontcorr,metrics.gefficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,18)
scatter(metrics.avecontcorr,metrics.defficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.avecontcorr,metrics.defficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel('AC v WD Correlation')

% MCR vs GT Metrics

subplot(4,5,4)
scatter(metrics.modalcontcorr,metrics.meanwdeg,50,metrics.group,'filled')
hold on
a = fit(metrics.modalcontcorr,metrics.meanwdeg,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
title('MC v WD Correlation')
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,9)
scatter(metrics.modalcontcorr,metrics.modularity,50,metrics.group,'filled')
hold on
a = fit(metrics.modalcontcorr,metrics.modularity,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,14)
scatter(metrics.modalcontcorr,metrics.gefficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.modalcontcorr,metrics.gefficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,19)
scatter(metrics.modalcontcorr,metrics.defficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.modalcontcorr,metrics.defficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel('MC v WD Correlation')

% AvM vs GT Metrics

subplot(4,5,5)
scatter(metrics.avmcorr,metrics.meanwdeg,50,metrics.group,'filled')
hold on
a = fit(metrics.avmcorr,metrics.meanwdeg,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
title('AC v MC Correlation')
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,10)
scatter(metrics.avmcorr,metrics.modularity,50,metrics.group,'filled')
hold on
a = fit(metrics.avmcorr,metrics.modularity,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,15)
scatter(metrics.avmcorr,metrics.gefficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.avmcorr,metrics.gefficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
xlabel([])
set(gca,'FontSize',10)

subplot(4,5,20)
scatter(metrics.avmcorr,metrics.defficiency,50,metrics.group,'filled')
hold on
a = fit(metrics.avmcorr,metrics.defficiency,'poly1');
b = plot(a,'predfunc')
legend('off');
b(1).Color = 'k';
b(2).Color = 'k';
b(3).Color = 'k';
set(gca,'FontSize',10)
xlabel('AC v MC Correlation')

%% Save

saveas(gcf,'FinalFigs/SupplementaryFigureS2.png')

