%% Load Connectomes

cd /Controllability;  % input location of data structures
patients = load('Patients/controllabilities.mat');
controls = load('Controls/controllabilities.mat'); 
vns = load('VNS/controllabilities.mat');
load('labels.mat')

%% Get ranks of wdeg, avecont and modalcont for each patient & control

for a = 1:length(patients.connectome)
    
    % wdegrank
    
    [~,r] = sort(patients.connectome(a).wdeg);
    rankV(r) = 1:numel(patients.connectome(a).wdeg);
    patients.connectome(a).wdegrank = (rankV)';
    
    % avecontrank
    
    [~,r] = sort(patients.connectome(a).avecont);
    rankV(r) = 1:numel(patients.connectome(a).avecont);
    patients.connectome(a).avecontrank = (rankV)';
    
    % modalcontrank
   
    [~,r] = sort(patients.connectome(a).modalcont);
    rankV(r) = 1:numel(patients.connectome(a).modalcont);
    patients.connectome(a).modalcontrank = (rankV)';
    
    % partcoefrank
   
%     [~,r] = sort(patients.connectome(a).partcoef);
%     rankV(r) = 1:numel(patients.connectome(a).partcoef);
%     patients.connectome(a).partcoefrank = (rankV)';
    
end

for a = 1:length(controls.connectome)
    
    % wdegrank
    
    [~,r] = sort(controls.connectome(a).wdeg);
    rankV(r) = 1:numel(controls.connectome(a).wdeg);
    controls.connectome(a).wdegrank = (rankV)';
    
    % avecontrank
    
    [~,r] = sort(controls.connectome(a).avecont);
    rankV(r) = 1:numel(controls.connectome(a).avecont);
    controls.connectome(a).avecontrank = (rankV)';
    
    % modalcontrank
   
    [~,r] = sort(controls.connectome(a).modalcont);
    rankV(r) = 1:numel(controls.connectome(a).modalcont);
    controls.connectome(a).modalcontrank = (rankV)';
    
    % partcoefrank
   
%     [~,r] = sort(controls.connectome(a).partcoef);
%     rankV(r) = 1:numel(controls.connectome(a).partcoef);
%     controls.connectome(a).partcoefrank = (rankV)';
%     
end

for a = 1:length(vns.connectome)
    
    % wdegrank
    
    [~,r] = sort(vns.connectome(a).wdeg);
    rankV(r) = 1:numel(vns.connectome(a).wdeg);
    vns.connectome(a).wdegrank = (rankV)';
    
    % avecontrank
    
    [~,r] = sort(vns.connectome(a).avecont);
    rankV(r) = 1:numel(vns.connectome(a).avecont);
    vns.connectome(a).avecontrank = (rankV)';
    
    % modalcontrank
   
    [~,r] = sort(vns.connectome(a).modalcont);
    rankV(r) = 1:numel(vns.connectome(a).modalcont);
    vns.connectome(a).modalcontrank = (rankV)';
    
    % partcoefrank
    
%     [~,r] = sort(vns.connectome(a).partcoef);
%     rankV(r) = 1:numel(vns.connectome(a).partcoef);
%     vns.connectome(a).partcoefrank = (rankV)';
    
end

%% Calculate correlation coefficients for each patient and control

for a = 1:length(controls.connectome)
    
    avecorr = corrcoef(controls.connectome(a).wdegrank,controls.connectome(a).avecontrank);
    modalcorr = corrcoef(controls.connectome(a).wdegrank,controls.connectome(a).modalcontrank);
    avmcorr = corrcoef(controls.connectome(a).avecontrank,controls.connectome(a).modalcontrank);

    controls.connectome(a).avecontcorr = avecorr(2);
    controls.connectome(a).modalcontcorr = modalcorr(2);
    controls.connectome(a).avmcorr = avmcorr(2);
    
    clear avecorr modalcorr avmcorr
end
for b = 1:length(patients.connectome)
   
    avecorr = corrcoef(patients.connectome(b).wdegrank,patients.connectome(b).avecontrank);
    modalcorr = corrcoef(patients.connectome(b).wdegrank,patients.connectome(b).modalcontrank);
    avmcorr = corrcoef(patients.connectome(b).avecontrank,patients.connectome(a).modalcontrank);


     patients.connectome(b).avecontcorr = avecorr(2);
     patients.connectome(b).modalcontcorr = modalcorr(2);
     patients.connectome(b).avmcorr = avmcorr(2);
     
     clear avecorr modalcorr avmcorr
end

for b = 1:length(vns.connectome)
   
    avecorr = corrcoef(vns.connectome(b).wdegrank,vns.connectome(b).avecontrank);
    modalcorr = corrcoef(vns.connectome(b).wdegrank,vns.connectome(b).modalcontrank);
    avmcorr = corrcoef(vns.connectome(b).avecontrank,vns.connectome(b).modalcontrank);

     vns.connectome(b).avecontcorr = avecorr(2);
     vns.connectome(b).modalcontcorr = modalcorr(2);
     vns.connectome(b).avmcorr = avmcorr(2);
     
     
     clear avecorr modalcorr avmcorr
end

%% Calculate rank means and SD for controls

controlwdegrankmean = mean([controls.connectome.wdegrank],2);
controlavecontrankmean = mean([controls.connectome.avecontrank],2);
controlmodalcontrankmean = mean([controls.connectome.modalcontrank],2);
controlwdegrankstd = std([controls.connectome.wdegrank],0,2);
controlavecontrankstd = std([controls.connectome.avecontrank],0,2);
controlmodalcontrankstd = std([controls.connectome.modalcontrank],0,2);

%% Z score ranks for each control & patient & VNS

for a = 1:length(controls.connectome)
    for b = 1:253
        zwdeg(b) = (controls.connectome(a).wdegrank(b) - controlwdegrankmean(b))/controlwdegrankstd(b);
        zavecont(b) = (controls.connectome(a).avecontrank(b) - controlavecontrankmean(b))/controlavecontrankstd(b);
        zmodalcont(b) = (controls.connectome(a).modalcontrank(b) - controlmodalcontrankmean(b))/controlmodalcontrankstd(b); 
    end
    
    controls.connectome(a).zwdeg = zwdeg';
    controls.connectome(a).zavecont = zavecont';
    controls.connectome(a).zmodalcont = zmodalcont';
    
    clear zwdeg zavecont zmodalcont;
    
end

for a = 1:length(patients.connectome)
    for b = 1:253
        zwdeg(b) = (patients.connectome(a).wdegrank(b) - controlwdegrankmean(b))/controlwdegrankstd(b);
        zavecont(b) = (patients.connectome(a).avecontrank(b) - controlavecontrankmean(b))/controlavecontrankstd(b);
        zmodalcont(b) = (patients.connectome(a).modalcontrank(b) - controlmodalcontrankmean(b))/controlmodalcontrankstd(b); 
    end
    
    patients.connectome(a).zwdeg = zwdeg';
    patients.connectome(a).zavecont = zavecont';
    patients.connectome(a).zmodalcont = zmodalcont';
    
    clear zwdeg zavecont zmodalcont;
    
end

for a = 1:length(vns.connectome)
    for b = 1:253
        zwdeg(b) = (vns.connectome(a).wdegrank(b) - controlwdegrankmean(b))/controlwdegrankstd(b);
        zavecont(b) = (vns.connectome(a).avecontrank(b) - controlavecontrankmean(b))/controlavecontrankstd(b);
        zmodalcont(b) = (vns.connectome(a).modalcontrank(b) - controlmodalcontrankmean(b))/controlmodalcontrankstd(b); 
    end
    
    vns.connectome(a).zwdeg = zwdeg';
    vns.connectome(a).zavecont = zavecont';
    vns.connectome(a).zmodalcont = zmodalcont';
    
    clear zwdeg zavecont zmodalcont;
    
end

%% Find parcels that have a |Z-score| > threshold

threshold = 3.1; % 
parcels = (1:253)';

% Patients 

grouplabels.patients.wdeghigh = label(mean([patients.connectome.zwdeg],2)>threshold);
grouplabels.patients.wdeglow = label(mean([patients.connectome.zwdeg],2)<-threshold);

grouplabels.patients.wdeghighRSN = RSN_assignment(mean([patients.connectome.zwdeg],2)>threshold);
grouplabels.patients.wdeglowRSN = RSN_assignment(mean([patients.connectome.zwdeg],2)<-threshold);

grouplabels.patients.aveconthigh = label(mean([patients.connectome.zavecont],2)>threshold);
grouplabels.patients.avecontlow = label(mean([patients.connectome.zavecont],2)<-threshold);

grouplabels.patients.aveconthighRSN = RSN_assignment(mean([patients.connectome.zavecont],2)>threshold);
grouplabels.patients.avecontlowRSN = RSN_assignment(mean([patients.connectome.zavecont],2)<-threshold);

grouplabels.patients.modalconthigh = label(mean([patients.connectome.zmodalcont],2)>threshold);
grouplabels.patients.modalcontlow = label(mean([patients.connectome.zmodalcont],2)<-threshold);

grouplabels.patients.modalconthighRSN = RSN_assignment(mean([patients.connectome.zmodalcont],2)>threshold);
grouplabels.patients.modalcontlowRSN = RSN_assignment(mean([patients.connectome.zmodalcont],2)<-threshold);

% VNS 

grouplabels.vns.wdeghigh = label(mean([vns.connectome.zwdeg],2)>threshold);
grouplabels.vns.wdeglow = label(mean([vns.connectome.zwdeg],2)<-threshold);

grouplabels.vns.aveconthigh = label(mean([vns.connectome.zavecont],2)>threshold);
grouplabels.vns.avecontlow = label(mean([vns.connectome.zavecont],2)<-threshold);

grouplabels.vns.modalconthigh = label(mean([vns.connectome.zmodalcont],2)>threshold);
grouplabels.vns.modalcontlow = label(mean([vns.connectome.zmodalcont],2)<-threshold);

%VNS RSN

grouplabels.vns.wdeghighRSN = RSN_assignment(mean([vns.connectome.zwdeg],2)>threshold);
grouplabels.vns.wdeglowRSN = RSN_assignment(mean([vns.connectome.zwdeg],2)<-threshold);

grouplabels.vns.aveconthighRSN = RSN_assignment(mean([vns.connectome.zavecont],2)>threshold);
grouplabels.vns.avecontlowRSN = RSN_assignment(mean([vns.connectome.zavecont],2)<-threshold);

grouplabels.vns.modalconthighRSN = RSN_assignment(mean([vns.connectome.zmodalcont],2)>threshold);
grouplabels.vns.modalcontlowRSN = RSN_assignment(mean([vns.connectome.zmodalcont],2)<-threshold);

%% Calculate mean Z-scores for resected parcels & ROB vs controls

for a = 1:length(patients.connectome)
    
% Find mean Z values for whole brain

patients.connectome(a).resection.zwdegwb = mean(patients.connectome(a).zwdeg);
patients.connectome(a).resection.zavecontwb = mean(patients.connectome(a).zavecont);
patients.connectome(a).resection.zmodalcontwb = mean(patients.connectome(a).zmodalcont);


% Find mean Z value for the resected parcels in that patient

patients.connectome(a).resection.zwdeg = mean(patients.connectome(a).zwdeg(patients.connectome(a).resectedparcels));
patients.connectome(a).resection.zavecont = mean(patients.connectome(a).zavecont(patients.connectome(a).resectedparcels));
patients.connectome(a).resection.zmodalcont = mean(patients.connectome(a).zmodalcont(patients.connectome(a).resectedparcels));

% Find mean Z value for rest of brain

patients.connectome(a).resection.zwdegrob = mean(patients.connectome(a).zwdeg(not(ismember(parcels,patients.connectome(a).resectedparcels))));
patients.connectome(a).resection.zavecontrob = mean(patients.connectome(a).zavecont(not(ismember(parcels,patients.connectome(a).resectedparcels))));
patients.connectome(a).resection.zmodalcontrob = mean(patients.connectome(a).zmodalcont(not(ismember(parcels,patients.connectome(a).resectedparcels))));

end

%% Create 1000 virtual resections for each patient of the same parcel size, include cortical/Hip/amy parcels only

cortexonly = setdiff(1:253, [109:119 122:123 235:245 248:253]); % Only allows resection of cortex/hipp/amyg

for a = 1:length(patients.connectome)
    for b = 1:1000
        patients.connectome(a).resection.virtual(b).parcels = cortexonly(randperm(numel(cortexonly),length(patients.connectome(a).resectedparcels)));
        patients.connectome(a).resection.virtual(b).zwdeg = mean(patients.connectome(a).zwdeg(patients.connectome(a).resection.virtual(b).parcels));
        patients.connectome(a).resection.virtual(b).zavecont = mean(patients.connectome(a).zavecont(patients.connectome(a).resection.virtual(b).parcels));
        patients.connectome(a).resection.virtual(b).zmodalcont = mean(patients.connectome(a).zmodalcont(patients.connectome(a).resection.virtual(b).parcels));

    end
    
    patients.connectome(a).resection.zofzwdeg = (patients.connectome(a).resection.zwdeg - mean([patients.connectome(a).resection.virtual.zwdeg]))/std([patients.connectome(a).resection.virtual.zwdeg]);
    patients.connectome(a).resection.zofzavecont = (patients.connectome(a).resection.zavecont - mean([patients.connectome(a).resection.virtual.zavecont]))/std([patients.connectome(a).resection.virtual.zavecont]);
    patients.connectome(a).resection.zofzmodalcont = (patients.connectome(a).resection.zmodalcont - mean([patients.connectome(a).resection.virtual.zmodalcont]))/std([patients.connectome(a).resection.virtual.zmodalcont]);
    
end
