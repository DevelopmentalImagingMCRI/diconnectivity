function M = mrtrix_create_tracks_ipsi_rois_crossing(S, SeedIMG, StartEndRegions, ShortLabels)

% marks the tracks that connect ipsilateral regions that go into the other 
% hemisphere's white matter or that connect interhemispheric regions that
% only go into one hemisphere's white matter as 1

%M = false(size(S));

LeftWMLabel = 2;
RightWMLabel = 41;

ShortLabelsHemis = cellfun(@(x) (x(1)), ShortLabels);
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
IMGLabels = SeedIMG(IDX);
clear IDX;

LeftWMLabels = (IMGLabels == LeftWMLabel);
RightWMLabels = (IMGLabels == RightWMLabel);

LeftWMLabels = mat2cell_vec(LeftWMLabels, TracksSZ);
RightWMLabels = mat2cell_vec(RightWMLabels, TracksSZ);

StreamlineHemis = ShortLabelsHemis(StartEndRegions);

M = false(length(S), 1);

LeftIpsiMask = all(StreamlineHemis == 'l', 2);
RightIpsiMask = all(StreamlineHemis == 'r', 2);

if(any(LeftIpsiMask))
	M(LeftIpsiMask) = cellfun(@any, RightWMLabels(LeftIpsiMask));
end
if(any(RightIpsiMask))
	M(RightIpsiMask) = cellfun(@any, LeftWMLabels(RightIpsiMask));
end
% OldM = M;
% 
% IMGLabels = mat2cell_vec(IMGLabels, TracksSZ);
% M = false(length(IMGLabels), 1);
% 
% for z = 1:length(IMGLabels)
% 	
% 	switch(StreamlineHemis(z, :))
% 		case 'll'
% 			% I am looking for label 2 but no label 41
% 			M(z) = any(IMGLabels{z} == RightWMLabel);
% 		case 'rr'
% 			M(z) = any(IMGLabels{z} == LeftWMLabel);
% 	end
% end
% 
% isequal(OldM, M)