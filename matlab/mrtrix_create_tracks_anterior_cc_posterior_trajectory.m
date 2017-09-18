function [M] = mrtrix_create_tracks_anterior_cc_posterior_trajectory(S, SeedIMG)

%M = false(size(S));

%C = 517;
% 253 is the central CC label
I = find(SeedIMG == 253);
[ID, ~, ~] = ind2sub(size(SeedIMG), I);
CentralCCCentroid = (min(ID) + max(ID)) / 2; %mean(ID)

%S{C}(:, 2)'
SZSlice = int32(size(SeedIMG, 1) * size(SeedIMG, 2));
SZRows = int32(size(SeedIMG, 1));

T = cat(1, S{:});
TracksSZ = cellfun('size', S, 1);

PastCCCentroid = (T(:, 2) > CentralCCCentroid);
PastCCCentroid = mat2cell_vec(PastCCCentroid, TracksSZ);
PastCCCentroid = cellfun(@any, PastCCCentroid);

T = int32(round(T));
IDX = SZSlice * (T(:, 3) - 1) + SZRows * (T(:, 1) - 1) + T(:, 2);


SeedValues = SeedIMG(IDX);
SeedValuesCC = (SeedValues == 255);
clear SeedValues;
SeedValuesCC = mat2cell_vec(SeedValuesCC, TracksSZ);

SeedValuesCC = cellfun(@any, SeedValuesCC);

M = SeedValuesCC & PastCCCentroid;
