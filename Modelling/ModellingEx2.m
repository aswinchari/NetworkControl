%% Load Connectomes

cd /Controllability;  % input location of data structures
patients = load('Patients/controllabilities.mat');
controls = load('Controls/controllabilities.mat'); 
vns = load('VNS/controllabilities.mat');
load('labels.mat');

%% Find average and std matrices for each cohort

z = [controls.connectome.connectome];
y = reshape(z,253,253,16);

group(1).name = 'controls';
group(1).matrices = y;
group(1).average = mean(y,3);
group(1).std = std(y,0,3);

clear y z

z = [patients.connectome.connectome];
y = reshape(z,253,253,52);

group(2).name = 'patients';
group(2).matrices = y;
group(2).average = mean(y,3);
group(2).std = std(y,0,3);

clear y z

z = [vns.connectome.connectome];
y = reshape(z,253,253,27);

group(3).name = 'vns';
group(3).matrices = y;
group(3).average = mean(y,3);
group(3).std = std(y,0,3);

clear y z 

z = [];
for a = 1:length(patients.connectome)
    if patients.connectome(a).atl == 1
        z = horzcat(z,patients.connectome(a).connectome);
    end
end
y = reshape(z,253,253,14);

group(4).name = 'patients only ATL';
group(4).matrices = y;
group(4).average = mean(y,3);
group(4).std = std(y,0,3);

clear y z 

z = [];
for a = 1:length(patients.connectome)
    if patients.connectome(a).atl == 0
        z = horzcat(z,patients.connectome(a).connectome);
    end
end
y = reshape(z,253,253,38);

group(5).name = 'patients not ATL';
group(5).matrices = y;
group(5).average = mean(y,3);
group(5).std = std(y,0,3);

clear y z 

%% Plot some average matrices = raw numbers

q = (group(2).average-group(1).average);
r = (group(3).average-group(1).average);

subplot(6,4,[1 2 5 6])
imagesc(q)
%set(gca,'ColorScale','log')
title('Resective Surgery - Controls')
colorbar
set(gca,'FontSize',12);
caxis([-200 200])

subplot(6,4,[3 4 7 8])
imagesc(r)
%set(gca,'ColorScale','log')
title('VNS - Controls')
set(gca,'FontSize',12);
colorbar
caxis([-200 200])

x = 1:1000;


for a = 1:length(x)
    y1high(a) = sum(sum(q>x(a)));
    y1low(a) = sum(sum(q<-x(a)));
    y2high(a) = sum(sum(r>x(a)));
    y2low(a) = sum(sum(r<-x(a)));    
end 

subplot(6,4,[9 10 13 14])
hold on
plot(x,y1high,'LineWidth',2,'Color','#77AC30')
plot(x,y1low,'LineWidth',2,'Color','#A2142F')
title('Difference in Edge Weights')
set(gca,'FontSize',12);
xlim([0 200])
ylim([0 10000])
xlabel('Threshold of Difference')
ylabel('Number of Edges')
legend('Resective Surgery > Controls','Controls > Resective Surgery')

subplot(6,4,[11 12 15 16])
hold on
plot(x,y2high,'LineWidth',2,'Color','#77AC30')
plot(x,y2low,'LineWidth',2,'Color','#A2142F')
title('Difference in Edge Weights')
set(gca,'FontSize',12);
xlim([0 200])
ylim([0 10000])
xlabel('Threshold of Difference')
ylabel('Number of Edges')
legend('VNS > Controls','Controls > VNS')

subplot(6,4,21)
imagesc(and(q>20,q<80))
title('Resective Surgery > Controls')
xticklabels([])
yticklabels([])

subplot(6,4,22)
imagesc(and(q<-20,q>-80))
title('Controls > Resective Surgery')
xticklabels([])
yticklabels([])

subplot(6,4,23)
imagesc(and(r>20,q<80))
title('VNS > Controls')
xticklabels([])
yticklabels([])

subplot(6,4,24)
imagesc(and(r<-20,q>-80))
title('Controls > VNS')
xticklabels([])
yticklabels([])

%% Save

saveas(gcf,'FinalFigs/groupdifferences.png')

%% Add weight (20-80 streamlines) just in the thalamocortical connections vs null models

for x = 1:length(controls.connectome)
    
    disp(strcat('ControlSubject:',string(x)))
    
    for d = 1:10
    
        nodesfromr = randi([109,123],1,12);
        nodestor = randi([1,108],1,60);
        nodesfroml = randi([235,249],1,12);
        nodestol = randi([124,234],1,60);
    
    z = controls.connectome(x).connectome;
    
    controls.connectome(x).model(211).connectome = z;
    
    for b = 1:length(nodesfromr)
        for c = 1:length(nodestor)
            
            aw = normrnd(50,15);
            
            z(nodesfromr(b),nodestor(c)) = z(nodesfromr(b),nodestor(c)) + aw;
            z(nodestor(c),nodesfromr(b)) = z(nodestor(c),nodesfromr(b)) + aw;
        end
    end
   
    for b = 1:length(nodesfroml)
        for c = 1:length(nodestol)
            
            aw = normrnd(50,15);
            
            z(nodesfroml(b),nodestol(c)) = z(nodesfroml(b),nodestol(c)) + aw;
            z(nodestol(c),nodesfroml(b)) = z(nodestol(c),nodesfroml(b)) + aw;
        end
    end
    
    z(z<0) = 0;
    controls.connectome(x).model(200+d).connectome = z;
    
    end
    
    disp('Making Matrices')
    
    for d = 1:200

            nodesfromr = randi([1,123],1,12);
            nodestor = randi([1,123],1,60);
            nodesfroml = randi([124,249],1,12);
            nodestol = randi([124,249],1,60);
    
    z = controls.connectome(x).connectome;
    
    for b = 1:length(nodesfromr)
        for c = 1:length(nodestor)
            
            aw = normrnd(50,15);
            
            z(nodesfromr(b),nodestor(c)) = z(nodesfromr(b),nodestor(c)) + aw;
            z(nodestor(c),nodesfromr(b)) = z(nodestor(c),nodesfromr(b)) + aw;
        end
    end
   
    for b = 1:length(nodesfroml)
        for c = 1:length(nodestol)
            
            aw = normrnd(50,15);
            
            z(nodesfroml(b),nodestol(c)) = z(nodesfroml(b),nodestol(c)) + aw;
            z(nodestol(c),nodesfroml(b)) = z(nodestol(c),nodesfroml(b)) + aw;
        end
    end
   
        z(z<0) = 0;
        controls.connectome(x).model(d).connectome = z;
        
    end
    
    disp('Calculating Metrics')
    
    for a = 1:length(controls.connectome(x).model)

    A = controls.connectome(x).model(a).connectome;
    NormA = A./(1+svds(A,1));     % Matrix normalization 
    [U, T] = schur(NormA,'real'); % Schur stability
    
    % Calculaute avecont 
    
    midMat = (U.^2)';
    v = diag(T);
    P = repmat(diag(1 - v*v'),1,size(NormA,1));
    controls.connectome(x).model(a).avecont = sum(midMat./P)';

    % Calculate modalcont 
    
    eigVals = diag(T);
    N = size(NormA,1);
    phi = zeros(N,1);
    for i = 1 : N
        phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
    end
    controls.connectome(x).model(a).modalcont = phi;
    
    % Weighted degree
    
    controls.connectome(x).model(a).wdeg = sum(A,2);
    
    % Clear parameters
    
    clear A NormA T U midMat v P eigVals N phi r rankV i
    
    % wdegrank
    
    [~,r] = sort(controls.connectome(x).model(a).wdeg);
    rankV(r) = 1:numel(controls.connectome(x).model(a).wdeg);
    controls.connectome(x).model(a).wdegrank = (rankV)';
    
    % avecontrank
    
    [~,r] = sort(controls.connectome(x).model(a).avecont);
    rankV(r) = 1:numel(controls.connectome(x).model(a).avecont);
    controls.connectome(x).model(a).avecontrank = (rankV)';
    
    % modalcontrank
   
    [~,r] = sort(controls.connectome(x).model(a).modalcont);
    rankV(r) = 1:numel(controls.connectome(x).model(a).modalcont);
    controls.connectome(x).model(a).modalcontrank = (rankV)';
   
    clear  r rankV
    
    % Graph theory metrics
    
    [module,modularity] = community_louvain(controls.connectome(x).model(a).connectome);
    gefficiency = efficiency_wei(controls.connectome(x).model(a).connectome);
    [defficiency,~] = diffusion_efficiency(controls.connectome(x).model(a).connectome);
   
    controls.connectome(x).model(a).module = module;
    controls.connectome(x).model(a).modularity = modularity;
    controls.connectome(x).model(a).gefficiency = gefficiency;
    controls.connectome(x).model(a).defficiency = defficiency;
   
    clear modularity gefficiency defficiency module
    
    AvMCorr = corrcoef(controls.connectome(x).model(a).avecontrank,controls.connectome(x).model(a).modalcontrank);
    AveContCorr = corrcoef(controls.connectome(x).model(a).wdegrank,controls.connectome(x).model(a).avecontrank);
    ModalContCorr = corrcoef(controls.connectome(x).model(a).wdegrank,controls.connectome(x).model(a).modalcontrank);
    
    controls.connectome(x).model(a).AvMCorr = AvMCorr(2);
    controls.connectome(x).model(a).AveContCorr = AveContCorr(2);
    controls.connectome(x).model(a).ModalContCorr = ModalContCorr(2);
    
    clear AvMCorr AveContCorr ModalContCorr
    
end
    
end

%% Plot that

    for a = 1:length(controls.connectome)
        
        for c = 1:211
        
      % Z-score for all models (200 null, 10 thalamic, baseline)
        
        controls.connectome(a).modelz(c).wdeg = (mean(controls.connectome(a).model(c).wdeg) - mean(mean([controls.connectome(a).model(1:200).wdeg])))/std(mean([controls.connectome(a).model(1:200).wdeg]));
        controls.connectome(a).modelz(c).avecont = (mean(controls.connectome(a).model(c).avecont) - mean(mean([controls.connectome(a).model(1:200).avecont])))/std(mean([controls.connectome(a).model(1:200).avecont]));
        controls.connectome(a).modelz(c).modalcont = (mean(controls.connectome(a).model(c).modalcont) - mean(mean([controls.connectome(a).model(1:200).modalcont])))/std(mean([controls.connectome(a).model(1:200).modalcont]));
        controls.connectome(a).modelz(c).modularity = (controls.connectome(a).model(c).modularity - mean([controls.connectome(a).model(1:200).modularity]))/std([controls.connectome(a).model(1:200).modularity]);
        controls.connectome(a).modelz(c).gefficiency = (controls.connectome(a).model(c).gefficiency - mean([controls.connectome(a).model(1:200).gefficiency]))/std([controls.connectome(a).model(1:200).gefficiency]);
        controls.connectome(a).modelz(c).defficiency = (controls.connectome(a).model(c).defficiency - mean([controls.connectome(a).model(1:200).defficiency]))/std([controls.connectome(a).model(1:200).defficiency]);
        controls.connectome(a).modelz(c).AveContCorr = (controls.connectome(a).model(c).AveContCorr - mean([controls.connectome(a).model(1:200).AveContCorr]))/std([controls.connectome(a).model(1:200).AveContCorr]);
        controls.connectome(a).modelz(c).ModalContCorr = (controls.connectome(a).model(c).ModalContCorr - mean([controls.connectome(a).model(1:200).ModalContCorr]))/std([controls.connectome(a).model(1:200).ModalContCorr]);
        controls.connectome(a).modelz(c).AvMCorr = (controls.connectome(a).model(c).AvMCorr - mean([controls.connectome(a).model(1:200).AvMCorr]))/std([controls.connectome(a).model(1:200).AvMCorr]);
        end
        
            
    end
    
 % Structure for all z scores   
    
    zscores = {};
    
for a = 1:16
    
    if a == 1
        zscores{1} = [controls.connectome(a).modelz(1:200).avecont]';
        zscores{2} = [controls.connectome(a).modelz(1:200).modalcont]';
        zscores{3} = [controls.connectome(a).modelz(1:200).AveContCorr]';
        zscores{4} = [controls.connectome(a).modelz(1:200).ModalContCorr]';
        zscores{5} = [controls.connectome(a).modelz(1:200).AvMCorr]';
        zscores{6} = [controls.connectome(a).modelz(1:200).modularity]';
        zscores{7} = [controls.connectome(a).modelz(1:200).gefficiency]';
        zscores{8} = [controls.connectome(a).modelz(1:200).defficiency]';
        
    else
        zscores{1} = vertcat(zscores{1}, [controls.connectome(a).modelz(1:200).avecont]');
        zscores{2} = vertcat(zscores{2}, [controls.connectome(a).modelz(1:200).modalcont]');
        zscores{3} = vertcat(zscores{3}, [controls.connectome(a).modelz(1:200).AveContCorr]');
        zscores{4} = vertcat(zscores{4}, [controls.connectome(a).modelz(1:200).ModalContCorr]');
        zscores{5} = vertcat(zscores{5}, [controls.connectome(a).modelz(1:200).AvMCorr]');
        zscores{6} = vertcat(zscores{6}, [controls.connectome(a).modelz(1:200).modularity]');
        zscores{7} = vertcat(zscores{7}, [controls.connectome(a).modelz(1:200).gefficiency]');
        zscores{8} = vertcat(zscores{8}, [controls.connectome(a).modelz(1:200).defficiency]');
    end
    
end

% Plot    

jit = (rand(10,1)-0.5)/4; 
cols = cbrewer('qual', 'Set1', 8); 

hold on
violin(zscores,'mc',[],'medc',[],'facecolor',cols,'facealpha',0.3);
% thalamic addition
for a = 1:length(controls.connectome)
    scatter(ones(10,1)*1+jit,[controls.connectome(a).modelz(201:210).avecont],20,cols(1,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*2+jit,[controls.connectome(a).modelz(201:210).modalcont],20,cols(2,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*3+jit,[controls.connectome(a).modelz(201:210).AveContCorr],20,cols(3,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*4+jit,[controls.connectome(a).modelz(201:210).ModalContCorr],20,cols(4,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*5+jit,[controls.connectome(a).modelz(201:210).AvMCorr],20,cols(5,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*6+jit,[controls.connectome(a).modelz(201:210).modularity],20,cols(6,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*7+jit,[controls.connectome(a).modelz(201:210).gefficiency],20,cols(7,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    scatter(ones(10,1)*8+jit,[controls.connectome(a).modelz(201:210).defficiency],20,cols(8,:),'filled','MarkerEdgeColor','k','MarkerFaceAlpha','flat','AlphaData',repelem(0.6,10))
    
    scatter(1-0.5,[controls.connectome(a).modelz(211).avecont],500,cols(1,:),'filled','p','MarkerEdgeColor','k')
    scatter(2-0.5,[controls.connectome(a).modelz(211).modalcont],500,cols(2,:),'filled','p','MarkerEdgeColor','k')
    scatter(3-0.5,[controls.connectome(a).modelz(211).AveContCorr],500,cols(3,:),'filled','p','MarkerEdgeColor','k')
    scatter(4-0.5,[controls.connectome(a).modelz(211).ModalContCorr],500,cols(4,:),'filled','p','MarkerEdgeColor','k')
    scatter(5-0.5,[controls.connectome(a).modelz(211).AvMCorr],500,cols(5,:),'filled','p','MarkerEdgeColor','k')
    scatter(6-0.5,[controls.connectome(a).modelz(211).modularity],500,cols(6,:),'filled','p','MarkerEdgeColor','k')
    scatter(7-0.5,[controls.connectome(a).modelz(211).gefficiency],500,cols(7,:),'filled','p','MarkerEdgeColor','k')
    scatter(8-0.5,[controls.connectome(a).modelz(211).defficiency],500,cols(8,:),'filled','p','MarkerEdgeColor','k')
end
ylim([-6 10])
yline(0);
xlim([0 8.5])
ylabel('Z-score')
xticklabels([])
xticks([])
box off

set(gca,'FontSize',15)
set(gca,'Xcolor','none')


%% Save

saveas(gcf,'FinalFigs/thalamocortical.png')

%% Did the metrics go up or down?

% set up constructs

% null models
    baselinen.avecont = [];
    baselinen.modalcont = [];
    baselinen.AveContCorr = [];
    baselinen.ModalContCorr = [];
    baselinen.AvMCorr = [];
    baselinen.modularity = [];
    baselinen.gefficiency = [];
    baselinen.defficiency = [];

for a = 1:length(controls.connectome)
    baselinen.avecont = [baselinen.avecont [controls.connectome(a).modelz(1:200).avecont]];
    baselinen.modalcont = [baselinen.modalcont [controls.connectome(a).modelz(1:200).modalcont]];
    baselinen.AveContCorr = [baselinen.AveContCorr [controls.connectome(a).modelz(1:200).AveContCorr]];
    baselinen.ModalContCorr = [baselinen.ModalContCorr [controls.connectome(a).modelz(1:200).ModalContCorr]];
    baselinen.AvMCorr = [baselinen.AvMCorr [controls.connectome(a).modelz(1:200).AvMCorr]];
    baselinen.modularity = [baselinen.modularity [controls.connectome(a).modelz(1:200).modularity]];
    baselinen.gefficiency = [baselinen.gefficiency [controls.connectome(a).modelz(1:200).gefficiency]];
    baselinen.defficiency = [baselinen.defficiency [controls.connectome(a).modelz(1:200).defficiency]];
end

% baseline

for a = 1:length(controls.connectome)
    baseline(a).avecont = controls.connectome(a).modelz(211).avecont;
    baseline(a).modalcont = controls.connectome(a).modelz(211).modalcont;
    baseline(a).AveContCorr = controls.connectome(a).modelz(211).AveContCorr;
    baseline(a).ModalContCorr = controls.connectome(a).modelz(211).ModalContCorr;
    baseline(a).AvMCorr = controls.connectome(a).modelz(211).AvMCorr;
    baseline(a).modularity = controls.connectome(a).modelz(211).modularity;
    baseline(a).gefficiency = controls.connectome(a).modelz(211).gefficiency;
    baseline(a).defficiency = controls.connectome(a).modelz(211).defficiency;
end

% thalamic models

    baselinet.avecont = [];
    baselinet.modalcont = [];
    baselinet.AveContCorr = [];
    baselinet.ModalContCorr = [];
    baselinet.AvMCorr = [];
    baselinet.modularity = [];
    baselinet.gefficiency = [];
    baselinet.defficiency = [];

for a = 1:length(controls.connectome)
    baselinet.avecont = [baselinet.avecont [controls.connectome(a).modelz(201:210).avecont]];
    baselinet.modalcont = [baselinet.modalcont [controls.connectome(a).modelz(201:210).modalcont]];
    baselinet.AveContCorr = [baselinet.AveContCorr [controls.connectome(a).modelz(201:210).AveContCorr]];
    baselinet.ModalContCorr = [baselinet.ModalContCorr [controls.connectome(a).modelz(201:210).ModalContCorr]];
    baselinet.AvMCorr = [baselinet.AvMCorr [controls.connectome(a).modelz(201:210).AvMCorr]];
    baselinet.modularity = [baselinet.modularity [controls.connectome(a).modelz(201:210).modularity]];
    baselinet.gefficiency = [baselinet.gefficiency [controls.connectome(a).modelz(201:210).gefficiency]];
    baselinet.defficiency = [baselinet.defficiency [controls.connectome(a).modelz(201:210).defficiency]];
end

% ttest for controls vs thalamic

fields = fieldnames(baseline);

for a = 1:length(fields)
    test(a).name = fields(a);
end

[~,test(1).p] = ttest2([baseline.avecont],[baselinet.avecont]);
[~,test(2).p] = ttest2([baseline.modalcont],[baselinet.modalcont]);
[~,test(3).p] = ttest2([baseline.AveContCorr],[baselinet.AveContCorr]);
[~,test(4).p] = ttest2([baseline.ModalContCorr],[baselinet.ModalContCorr]);
[~,test(5).p] = ttest2([baseline.AvMCorr],[baselinet.AvMCorr]);
[~,test(6).p] = ttest2([baseline.modularity],[baselinet.modularity]);
[~,test(7).p] = ttest2([baseline.gefficiency],[baselinet.gefficiency]);
[~,test(8).p] = ttest2([baseline.defficiency],[baselinet.defficiency]);

[correction,~,~] = fdr_BH([test.p],0.05);

for a = 1:8
    test(a).corrp = correction(a);
end

% two sample ttests for controls vs null

fields = fieldnames(baselinet);

for a = 1:length(fields)
    testt(a).name = fields(a);
end

[~,testt(1).p] = ttest2([baseline.avecont],[baselinen.avecont]);
[~,testt(2).p] = ttest2([baseline.modalcont],[baselinen.modalcont]);
[~,testt(3).p] = ttest2([baseline.AveContCorr],[baselinen.AveContCorr]);
[~,testt(4).p] = ttest2([baseline.ModalContCorr],[baselinen.ModalContCorr]);
[~,testt(5).p] = ttest2([baseline.AvMCorr],[baselinen.AvMCorr]);
[~,testt(6).p] = ttest2([baseline.modularity],[baselinen.modularity]);
[~,testt(7).p] = ttest2([baseline.gefficiency],[baselinen.gefficiency]);
[~,testt(8).p] = ttest2([baseline.defficiency],[baselinen.defficiency]);


[correction,~,~] = fdr_BH([testt.p],0.05);

for a = 1:8
    testt(a).corrp = correction(a);
end

