%% Load Connectomes

cd /Controllability;  % input location of data structures
controls = load('Controls/controllabilities.mat'); 

%% Stipulate how many edges to add

SynthGraphs = 200;                                  % number of synthetic graphs
AddEdges = 1000;                                    % number of edges to add (1000 ~= 3%)
pd = makedist('Weibull','a',0.717,'b',0.325);       % distribution of edge weights

%% Manipulate the base model

for x = 1:length(controls.connectome)
    
    disp(strcat('Subject',num2str(x)))
    
  for a = 1:SynthGraphs   
   
    avg = controls.connectome(x).connectome;
    [ro co] = find(avg==0); % find zero elements
      
    for b = 1:AddEdges      % Adding Edges
        
        c = randi(length(ro));
        avg(ro(c),co(c)) = random(pd);
        avg(co(c),ro(c)) = avg(ro(c),co(c)); 
        
    end
    
    avg = avg - diag(diag(avg)); % ensure diagonals are zero
    
    controls.connectome(x).model(a).connectome = avg;
    controls.connectome(x).model(a).name = a;
    controls.connectome(x).model(a).colour = 2;

clear avg

end
%% Add the real connectome

controls.connectome(x).model(SynthGraphs+1).connectome = controls.connectome(x).connectome;
controls.connectome(x).model(SynthGraphs+1).name = 'base';
controls.connectome(x).model(SynthGraphs+1).colour = 1;

%% Calculate metrics

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

%% Plot that construct

% Z score

    for a = 1:length(controls.connectome)
        controls.connectome(a).modelwdeg = (mean(controls.connectome(a).model(201).wdeg) - mean(mean([controls.connectome(a).model.wdeg])))/std(mean([controls.connectome(a).model.wdeg]));
        controls.connectome(a).modelavecont = (mean(controls.connectome(a).model(201).avecont) - mean(mean([controls.connectome(a).model.avecont])))/std(mean([controls.connectome(a).model.avecont]));
        controls.connectome(a).modelmodalcont = (mean(controls.connectome(a).model(201).modalcont) - mean(mean([controls.connectome(a).model.modalcont])))/std(mean([controls.connectome(a).model.modalcont]));
        controls.connectome(a).modelmodularity = (controls.connectome(a).model(201).modularity - mean([controls.connectome(a).model.modularity]))/std([controls.connectome(a).model.modularity]);
        controls.connectome(a).modelgefficiency = (controls.connectome(a).model(201).gefficiency - mean([controls.connectome(a).model.gefficiency]))/std([controls.connectome(a).model.gefficiency]);
        controls.connectome(a).modeldefficiency = (controls.connectome(a).model(201).defficiency - mean([controls.connectome(a).model.defficiency]))/std([controls.connectome(a).model.defficiency]);
        controls.connectome(a).modelAveContCorr = (controls.connectome(a).model(201).AveContCorr - mean([controls.connectome(a).model.AveContCorr]))/std([controls.connectome(a).model.AveContCorr]);
        controls.connectome(a).modelModalContCorr = (controls.connectome(a).model(201).ModalContCorr - mean([controls.connectome(a).model.ModalContCorr]))/std([controls.connectome(a).model.ModalContCorr]);
        controls.connectome(a).modelAvMCorr = (controls.connectome(a).model(201).AvMCorr - mean([controls.connectome(a).model.AvMCorr]))/std([controls.connectome(a).model.AvMCorr]);
    end
    
% Plot    

jit = (rand(length(controls.connectome),1)-0.5)/4; 
cols = cbrewer('qual', 'Set1', 10); 

scatter(ones(length(controls.connectome),1)+jit,[controls.connectome.modelwdeg],100,cols(10,:),'filled','MarkerEdgeColor','k')
hold on
scatter(ones(length(controls.connectome),1)*2+jit,[controls.connectome.modelavecont],100,cols(1,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*3+jit,[controls.connectome.modelmodalcont],100,cols(2,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*4+jit,[controls.connectome.modelAveContCorr],100,cols(3,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*5+jit,[controls.connectome.modelModalContCorr],100,cols(4,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*6+jit,[controls.connectome.modelAvMCorr],100,cols(5,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*7+jit,[controls.connectome.modelmodularity],100,cols(6,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*8+jit,[controls.connectome.modelgefficiency],100,cols(8,:),'filled','MarkerEdgeColor','k')
scatter(ones(length(controls.connectome),1)*9+jit,[controls.connectome.modeldefficiency],100,cols(9,:),'filled','MarkerEdgeColor','k')
ylim([-10 6])
yline(0);
xlim([0.5 9.5])
ylabel('Z-score')
xticklabels([])
xticks([])

set(gca,'FontSize',15)


