%% Extracts the CSV connectomes and saves them as mat files
% Needsto be repeated for each group of subjects

listofconnectomes = dir ('/rawdata/sub*/ses-preop/dti/mrtrix/newconnectome.csv')

for a = 1:length(listofconnectomes)
    
    file = listofconnectomes(a)
    load(strcat(file.folder,'/',file.name), '-ascii')
    save(strcat(file.folder(42:47),'_connectome.mat'),'connectome')
    clear file connectome
    
end
