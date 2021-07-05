%% Assembles the connectomes for each group of patients in a structure

%% Controls

filelist = dir('/Users/aswinchari/Desktop/NewConnectome/Controls/*.mat')

for a = 1:length(filelist)
    load(filelist(a).name);
    connectome(a).connectome = newconnectome;
    clear newconnectome
end

for b = 1:length(filelist)
   connectome(b).sub = filelist(b).name(1:6);
   connectome(b).connectome = connectome(b).connectome - (diag(diag(connectome(b).connectome)));
end

save('Controls/allconnectomes.mat','connectome');

%% Resective Surgery

filelist = dir('/Users/aswinchari/Desktop/NewConnectome/Patients/*.mat')

for a = 1:length(filelist)
    load(filelist(a).name);
    connectome(a).connectome = newconnectome;
    clear newconnectome
end

for b = 1:length(filelist)
   connectome(b).sub = filelist(b).name(1:6);
   connectome(b).connectome = connectome(b).connectome - (diag(diag(connectome(b).connectome)));
end

save('Patients/allconnectomes.mat','connectome');

%% VNS

filelist = dir('/Users/aswinchari/Desktop/NewConnectome/VNS/*.mat')

for a = 1:length(filelist)
    load(filelist(a).name);
    connectome(a).connectome = newconnectome;
    clear newconnectome
end

for b = 1:length(filelist)
   connectome(b).sub = filelist(b).name(1:6);
   connectome(b).connectome = connectome(b).connectome - (diag(diag(connectome(b).connectome)));
end

save('VNS/allconnectomes.mat','connectome');