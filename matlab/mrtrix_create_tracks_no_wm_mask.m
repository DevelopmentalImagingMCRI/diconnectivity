function M = mrtrix_create_tracks_no_wm_mask(S, SeedIMG)

% marks the tracks that do not go through white matter as 1

%M = false(size(S));

SZSlice = int32(size(SeedIMG, 1) * size(SeedIMG, 2));
SZRows = int32(size(SeedIMG, 1));

T = cat(1, S{:});
TracksSZ = cellfun('size', S, 1);

%SeedValues = interp3_linear_fast(double(SeedIMG), double(T(:, 1)), double(T(:, 2)), double(T(:, 3)), 'nearest');
%SeedValues = (SeedValues == 0);
%SeedValues = mat2cell_vec(SeedValues, TracksSZ);
T = int32(round(T));
IDX = SZSlice * (T(:, 3) - 1) + SZRows * (T(:, 1) - 1) + T(:, 2);
clear T;

%NotWM = (SeedIMG(IDX) == 0);
%NotWM = mat2cell_vec(NotWM, TracksSZ);
%M = cellfun(@all, WM);

% should be more efficient to check for any WM voxel rather than all NonWM
% voxels
WM = (SeedIMG(IDX) > 0);
WM = mat2cell_vec(WM, TracksSZ);
M = ~cellfun(@any, WM);
