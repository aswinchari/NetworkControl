% Convert the atlas to surface mz3 (need to do this once only)

nii_nii2atlas('FinalPaedBrainLausanne.nii'); % this saves the file as merge.mz3

% Use any output which is a 253x1 vector which assigns a value to each
% parcel

x = rand(253,1);

% Use AtlasStatMap to convert this to a map

AtlasStatMap('merge.mz3','x.mz3',[],x);

% Now open this file in surfice