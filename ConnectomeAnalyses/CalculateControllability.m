%% Calculates weighted degree, average and modal controllability for all subjects within the structure

%% Controls

load('Controls/allconnectomes.mat')

for a = 1:length(connectome)

    A = connectome(a).connectome;
    NormA = A./(1+svds(A,1));     % Matrix normalization 
    [U, T] = schur(NormA,'real'); % Schur stability
    
    % Calculaute avecont 
    
    midMat = (U.^2)';
    v = diag(T);
    P = repmat(diag(1 - v*v'),1,size(NormA,1));
    connectome(a).avecont = sum(midMat./P)';

    % Calculate modalcont 
    
    eigVals = diag(T);
    N = size(NormA,1);
    phi = zeros(N,1);
    for i = 1 : N
        phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
    end
    connectome(a).modalcont = phi;
    
    % Weighted degrees & rank of each node's streamlines 
    
    connectome(a).wdeg = sum(A,2);
    [~,r] = sort(connectome(a).wdeg);
    rankV(r) = 1:numel(connectome(a).wdeg);
    connectome(a).rank = (rankV)';
    
    % Clear parameters
    
    clear A NormA T U midMat v P eigVals N phi r rankV
    
end

% Calculate mean

for b = 1:length(connectome)
    
    connectome(b).meanavecont = mean(connectome(b).avecont);
    connectome(b).meanmodalcont = mean(connectome(b).modalcont);
end

save('Controls/controllabilities.mat','connectome')

%% Resective Surgery

load('Patients/allconnectomes.mat')

for a = 1:length(connectome)

    A = connectome(a).connectome;
    NormA = A./(1+svds(A,1));     % Matrix normalization 
    [U, T] = schur(NormA,'real'); % Schur stability
    
    % Calculaute avecont 
    
    midMat = (U.^2)';
    v = diag(T);
    P = repmat(diag(1 - v*v'),1,size(NormA,1));
    connectome(a).avecont = sum(midMat./P)';

    % Calculate modalcont 
    
    eigVals = diag(T);
    N = size(NormA,1);
    phi = zeros(N,1);
    for i = 1 : N
        phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
    end
    connectome(a).modalcont = phi;
    
    % Weighted degrees & rank of each node's streamlines 
    
    connectome(a).wdeg = sum(A,2);
    [~,r] = sort(connectome(a).wdeg);
    rankV(r) = 1:numel(connectome(a).wdeg);
    connectome(a).rank = (rankV)';
    
    % Clear parameters
    
    clear A NormA T U midMat v P eigVals N phi r rankV
    
end

% Calculate mean

for b = 1:length(connectome)
    
    connectome(b).meanavecont = mean(connectome(b).avecont);
    connectome(b).meanmodalcont = mean(connectome(b).modalcont);
end

save('Patients/controllabilities.mat','connectome')


%% VNS

load('VNS/allconnectomes.mat')

for a = 1:length(connectome)

    A = connectome(a).connectome;
    NormA = A./(1+svds(A,1));     % Matrix normalization 
    [U, T] = schur(NormA,'real'); % Schur stability
    
    % Calculaute avecont 
    
    midMat = (U.^2)';
    v = diag(T);
    P = repmat(diag(1 - v*v'),1,size(NormA,1));
    connectome(a).avecont = sum(midMat./P)';

    % Calculate modalcont 
    
    eigVals = diag(T);
    N = size(NormA,1);
    phi = zeros(N,1);
    for i = 1 : N
        phi(i) = (U(i,:).^2) * (1 - eigVals.^2);
    end
    connectome(a).modalcont = phi;
    
    % Weighted degrees & rank of each node's streamlines 
    
    connectome(a).wdeg = sum(A,2);
    [~,r] = sort(connectome(a).wdeg);
    rankV(r) = 1:numel(connectome(a).wdeg);
    connectome(a).rank = (rankV)';
    
    % Clear parameters
    
    clear A NormA T U midMat v P eigVals N phi r rankV
    
end

% Calculate mean

for b = 1:length(connectome)
    
    connectome(b).meanavecont = mean(connectome(b).avecont);
    connectome(b).meanmodalcont = mean(connectome(b).modalcont);
end

save('VNS/controllabilities.mat','connectome')
