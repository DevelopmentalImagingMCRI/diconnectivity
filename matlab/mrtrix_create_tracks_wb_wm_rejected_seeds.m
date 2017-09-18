function [varargout] = mrtrix_create_tracks_wb_wm_rejected_seeds(Subject, MrtrixMethod, SeedType, Curvature)

% same as mrtrix_create_tracks_wb_wm but saves the rejected seeds image
% mrtrix_create_tracks_wb_wm_rejected_seeds(Subject, MrtrixMethod, SeedType, Curvature)
%
% DESCRIPTION
%	Performs whole brain, white matter tractography. Creates connectivity
%	matrices multiple times and writes all of them out. Does not save the
%	track files. Performs tractography NumIters times.
%
% INPUT FILES
%	Subject/freesurfer_wm.nii.gz
%	
% OUTPUT FILES
%	Subject/['tracks_' SeedType '_' MrtrixMethod '_accepted_seeds_' num2str(Curvature) '.nii.gz']
%   Subject/['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']),
%		'FinalStartRegions', 
%		'FinalEndRegions',
%		'WeightedA',
%		'CountA', 
%		'SizeWeightedA', 
%		'LengthA',
%		'SeedSizes')
%	Subject/['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '.tck'
%	Subject/['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '.trk'

CortexLabels = load_freesurfer_cortex_labels;

WMMaskNII = load_nii(fullfile(Subject, 'freesurfer_wm'));
WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);

% this creates too many tracks for high res volumes
% create two tracks per voxel of 0.75 x 0.75 x 0.75
WMMaskNIIVoxVolume = prod(WMMaskNII.hdr.dime.pixdim(2:4));

DesiredVolume = 0.7 * 0.7 * 0.7;
TotalNumTracks = round(WMMaskNIIVoxVolume / DesiredVolume * sum(WMMaskIMG(:) > 0));
clear DesiredVolume;
TotalNumTracks = 50000;

NumTracksSoFar = 0;
FinalTracks = cell(1, TotalNumTracks);
FinalArcLengths = zeros(1, TotalNumTracks);
FinalStartRegions = zeros(1, TotalNumTracks);
FinalEndRegions = zeros(1, TotalNumTracks);
TracksHeader = [];

switch(SeedType)
	case 'aparc'
		CurLabels = CortexLabels.APARC;
		ExcludeFile = fullfile(Subject, ['exclude_' SeedType '.nii.gz']);
	case 'aparc.a2009s'
		CurLabels = CortexLabels.APARCOhNine;
		ExcludeFile = fullfile(Subject, 'exclude_aparcohnine.nii.gz');
	case 'aparc_ads'
		CurLabels = CortexLabels.APARCADS;
		ExcludeFile = fullfile(Subject, ['exclude_' SeedType '.nii.gz']);
	case 'wmparc'
		CurLabels = CortexLabels.WMPARC;
		ExcludeFile = fullfile(Subject, ['exclude_' SeedType '.nii.gz']);
	case 'subdivided1000'
		CurLabels = CortexLabels.SUBDIVIDED;
		ExcludeFile = fullfile(Subject, ['exclude_' SeedType '.nii.gz']);
	case 'wmsubdivided1000'
		CurLabels = CortexLabels.WMSUBDIVIDED;
		ExcludeFile = fullfile(Subject, ['exclude_' SeedType '.nii.gz']);
end

switch(SeedType)
	case {'aparc', 'aparc.a2009s', 'aparc_ads'}
		IncludeFile = fullfile(Subject, ['brain_non_wm_' SeedType '.nii.gz']);
	case {'wmparc', 'wmsubdivided1000'}
		IncludeFile = fullfile(Subject, ['wm_' SeedType '.nii.gz']);
end

if(exist(IncludeFile, 'file') ~= 2)
	error(['IncludeFile ' IncludeFile ' does not exist']);
end

OriginalSeedNII = load_nii(fullfile('..', 'register_freesurfer_to_mrtrix_highres_masks', Subject, CurLabels.SeedFile));
OriginalSeedIMG = flipdim(permute(OriginalSeedNII.img, [2 1 3]), 1);
SeedNII = load_nii(IncludeFile);
SeedIMG = flipdim(permute(SeedNII.img, [2 1 3]), 1);

disp('Saving seeds');
SeedOutputDir = fullfile(Subject, ['seeds_' SeedType]);
[~,~,~] = mkdir(SeedOutputDir);
% for z = 1:length(CurLabels.values)
% 	CurSeedIMG = uint8(SeedIMG == z);
% 	CurSeedNII = SeedNII;
% 	CurSeedNII.img = permute(flipdim(CurSeedIMG, 1), [2 1 3]);
% 	save_nii(CurSeedNII, fullfile(Subject, ['seeds_' SeedType], [CurLabels.shortlabels{z} '.nii.gz']));
% 	clear CurSeedIMG CurSeedNII;
% end
disp('Finished saving seeds');

SeedSizes = histc(SeedIMG(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;

AcceptedSeedsNII = WMMaskNII;
AcceptedSeedsIMG = zeros(size(SeedIMG), 'int32');

WeightedA = zeros(length(CurLabels.values));
LengthA = zeros(length(CurLabels.values));
CountA = zeros(length(CurLabels.values));

MrtrixDirectory = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9-save-accepted-print-rejected')];
BlockSize = 50000;

while(NumTracksSoFar < TotalNumTracks)
	MrtrixStep = 1;
	NumTracks = max(min((TotalNumTracks - NumTracksSoFar) * 2, 2000000), 1000);
	TrackFile = [tempname '.tck'];
	switch(lower(MrtrixMethod))
		case 'dt_stream'
	%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
	CommandString = sprintf('LD_LIBRARY_PATH=%s; %s DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -curvature %f -stop -minlength %f 2>/dev/null;', ...
		fullfile(MrtrixDirectory, 'lib'), ...
		fullfile(MrtrixDirectory, 'bin', 'streamtrack'), ...
		fullfile(Subject, 'dwi.nii.gz'), ...
		TrackFile, ...
		fullfile(Subject, 'freesurfer_wm.nii.gz'), ...
		ExcludeFile, ...
		IncludeFile, ...
		fullfile(Subject, 'mask.nii.gz'), NumTracks, NumTracks * 100, ...
		MrtrixStep, ...
		fullfile(Subject, 'grad.b'), ...
		Curvature, ...
		10);
		case {'sd_stream', 'sd_prob', 'sd_runge2', 'sd_runge4'}
	%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
	CommandString = sprintf('LD_LIBRARY_PATH=%s; %s %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -stop -curvature %f -minlength %f 2>&1;', ...
		fullfile(MrtrixDirectory, 'lib'), ...
		fullfile(MrtrixDirectory, 'bin', 'streamtrack'), ...
		upper(MrtrixMethod), ...
		fullfile(Subject, 'CSD.nii.gz'), ...
		TrackFile, ...
		fullfile(Subject, 'freesurfer_wm.nii.gz'), ...
		ExcludeFile, ...
		IncludeFile, ...
		fullfile(Subject, 'mask.nii.gz'), NumTracks, NumTracks * 100, ...
		MrtrixStep, ...
		fullfile(Subject, 'grad.b'), ...
		Curvature, ...
		10);
	end
	disp(CommandString);
	disp('Generating tracks');
	[~, result] = system(CommandString);
	disp('Finished generating tracks');
	MIFTracks = load_mif_tracks(TrackFile);
	%copyfile(TrackFile, 'tmp.tck');
	delete(TrackFile);
	CurTracks = MIFTracks.Tracks;
	TracksSZ = cellfun('size', CurTracks, 1);
	CurNumTracks = numel(CurTracks);
	ArcLengths = TracksSZ * MrtrixStep;

	if(isempty(TracksHeader))
		TracksHeader = MIFTracks;
		TracksHeader.Tracks = [];
	end
	MIFTracks.Tracks = [];

	tokens = textscan(result, '%s %f %f %f %d %d %d %d %d');
	
	MrtrixAcceptedTrackMask = strcmp(tokens{1}, 'accept:');
	
	AcceptedTokens = cell(size(tokens));
	RejectedTokens = cell(size(tokens));
	
	for z = 1:length(tokens)
		RejectedTokens{z} = tokens{z}(~MrtrixAcceptedTrackMask);
	end
	
	FirstCurNumTracksAcceptedIDX = find(MrtrixAcceptedTrackMask, CurNumTracks, 'first');
	
	for z = 1:length(AcceptedTokens)
		AcceptedTokens{z} = tokens{z}(FirstCurNumTracksAcceptedIDX);
	end
	clear FirstCurNumTracksAcceptedIDX;
	clear MIFTracks;

	% convert the tracks from world coordinates (mrtrix) to voxel
	% coordinates
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = (z:RightIDX);

		CurTracks(CurIDX) = tracks_world_to_img(CurTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
	end

	%AcceptedTrackMask = strcmp(tokens{1}, 'accept:');
	AcceptedTrackMask = true(CurNumTracks, 1);
	%disp([num2str(sum(AcceptedTrackMask)) ' tracks left after insideimg']);

	disp([num2str(length(AcceptedTokens{1})) ' mrtrix accepted tracks']);
	disp([num2str(length(RejectedTokens{1})) ' mrtrix rejected tracks']);

	StartRegion = AcceptedTokens{8};
	EndRegion = AcceptedTokens{9};

	AcceptedTrackMask(StartRegion == EndRegion) = 0;
	disp([num2str(sum(AcceptedTrackMask)) ' tracks left after same region cull']);
	NumBeforeNoWM = sum(AcceptedTrackMask);

	CurNumTracks = sum(AcceptedTrackMask);
	IDX = find(AcceptedTrackMask);
	RejectedTracksNoWM = false(CurNumTracks, 1);
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = z:RightIDX;
		RejectedTracksNoWM(CurIDX) = tracks_no_wm_mask(CurTracks(IDX(CurIDX)), WMMaskIMG);
	end
	clear IDX RightIDX CurIDX;

	AcceptedTrackMask(AcceptedTrackMask) = ~RejectedTracksNoWM;
	clear RejectedTracksNoWM;
	NumBeforeAnteriorPosterior = sum(AcceptedTrackMask);
	disp([num2str(sum(AcceptedTrackMask)) ' tracks left after no_wm']);

	CurNumTracks = sum(AcceptedTrackMask);
	IDX = find(AcceptedTrackMask);

	% reject posterior fibres coming out of the anterior Corpus
	% Callosum
	RejectedAnteriorCCPosterior = false(CurNumTracks, 1);
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = z:RightIDX;
		RejectedAnteriorCCPosterior(CurIDX) = tracks_anterior_cc_posterior_trajectory(CurTracks(IDX(CurIDX)), OriginalSeedIMG);
	end
	AcceptedTrackMask(AcceptedTrackMask) = ~RejectedAnteriorCCPosterior;

	clear RejectedAnteriorCCPosterior;
	disp([num2str(NumBeforeNoWM - NumBeforeAnteriorPosterior) ' culled by no wm']);
	disp([num2str(NumBeforeAnteriorPosterior - sum(AcceptedTrackMask)) ' culled by anterior posterior']);
	AcceptedTracks = CurTracks(AcceptedTrackMask);

	Seeds = cat(2, AcceptedTokens{2:4});
	AcceptedSeeds = Seeds(AcceptedTrackMask, :);

	clear Seeds;
	AcceptedCodes = AcceptedTokens{5};

	UniqueCodes = unique(AcceptedCodes);
	for z = 1:length(UniqueCodes)
		disp([num2str(sum(AcceptedCodes == UniqueCodes(z))) ' got code ' num2str(UniqueCodes(z))]);
	end

	%clear AcceptedTokens;

	disp([num2str(sum(~AcceptedTrackMask)) ' rejected tracks']);

	AcceptedArcLengths = ArcLengths(AcceptedTrackMask);
	AcceptedStartRegions = StartRegion(AcceptedTrackMask);
	AcceptedEndRegions = EndRegion(AcceptedTrackMask);
	clear RejectedTracks;

	AcceptedSeedsVox = tracks_world_to_img({AcceptedSeeds}, WMMaskNII, size(WMMaskIMG));
	AcceptedSeedsVox = AcceptedSeedsVox{1};

	I = round(AcceptedSeedsVox(:, 2)) + size(AcceptedSeedsIMG, 1) * round(AcceptedSeedsVox(:, 1)) + size(AcceptedSeedsIMG, 1) * size(AcceptedSeedsIMG, 2) * round(AcceptedSeedsVox(:, 3));
	AddToIMG = histc(I, 1:numel(AcceptedSeedsIMG));
	AcceptedSeedsIMG = AcceptedSeedsIMG + int32(reshape(AddToIMG, size(AcceptedSeedsIMG)));
	clear I AddToIMG AcceptedSeedsVox;

	% update the images of the rejected seeds
	% code 1 is hits exclusion region
	% code 2 is success (unidirectional tracking)
	% code 3 is ??? stop_when_included && !unidirectional && includedfirst && !doingseconddirection) return (3);
	% code 4 is success, however, if the streamline was rejected it means the streamline was too short
	% code 5 is the streamline was too long
	% code 7 is if the streamline exited the mask
	% code -2 is nan max_value
	% code -3 is low FOD amplitude
	% code -4 is curvature threshold exceeded
	
	RejectedSeeds = cat(2, RejectedTokens{2:4});
	RejectedCodes = RejectedTokens{5};
	
	RejectedSeedsVox = tracks_world_to_img({RejectedSeeds}, WMMaskNII, size(WMMaskIMG));
	RejectedSeedsVox = RejectedSeedsVox{1};
	
	keyboard;
	
	% fast way of adding the counts to the connectivity matrices
	% get the linear indices of the regions in the connectivity
	% matrices
	CurRegionsI = sub2ind(size(CountA), AcceptedStartRegions, AcceptedEndRegions);

	[CurRegionsISorted, CurRegionsISortedIIDX] = sort(CurRegionsI);

	% get the sizes and indices of 
	FirstIndices = [1; find(diff(CurRegionsISorted) > 0) + 1];
	CurSZ = diff([FirstIndices; length(CurRegionsI) + 1]);

	LengthAValues = mat2cell_vec(AcceptedArcLengths(CurRegionsISortedIIDX), int32(CurSZ));
	WeightedAValues = mat2cell_vec(1 ./ AcceptedArcLengths(CurRegionsISortedIIDX), int32(CurSZ));

	CurCountA = zeros(size(CountA));
	CurWeightedA = zeros(size(CountA));
	CurLengthA = zeros(size(CountA));

	CurCountA(CurRegionsISorted(FirstIndices)) = CurSZ;
	CurWeightedA(CurRegionsISorted(FirstIndices)) = cellfun(@sum, WeightedAValues);
	CurLengthA(CurRegionsISorted(FirstIndices)) = cellfun(@sum, LengthAValues);

	CurCountA = CurCountA + CurCountA';
	CurWeightedA = CurWeightedA + CurWeightedA';
	CurLengthA = CurLengthA + CurLengthA';

	CountA = CountA + CurCountA;
	WeightedA = WeightedA + CurWeightedA;
	LengthA = LengthA + CurLengthA;

	clear CurCountA CurWeightedA CurLengthA WeightedAValues LengthAValues CurSZ FirstIndices CurRegionsISorted CurRegionsISortedIIDX CurRegionsI;

	TracksToAdd = 1:min(length(AcceptedTracks), TotalNumTracks - NumTracksSoFar);%, length(CurTracks));

	FinalTracks(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedTracks(TracksToAdd);
	FinalArcLengths(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedArcLengths(TracksToAdd);

	FinalStartRegions(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedStartRegions(TracksToAdd);
	FinalEndRegions(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedEndRegions(TracksToAdd);
	NumTracksSoFar = NumTracksSoFar + length(TracksToAdd);
	clear CurTracks SeedValues RejectedTracks AcceptedTracks;
	disp(['Number of tracks generated so far: ' num2str(NumTracksSoFar) ' of ' num2str(TotalNumTracks)]);
end

save(fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalTracks', 'FinalArcLengths', 'FinalStartRegions', 'FinalEndRegions');

NumFinalTracks = numel(FinalTracks);
% block to save RAM
%FinalTracks = tracks_img_to_world(FinalTracks, WMMaskNII, size(WMMaskIMG));

for z = 1:BlockSize:NumFinalTracks
	RightIDX = min(z + BlockSize - 1, NumFinalTracks);
	CurIDX = (z:RightIDX);

	FinalTracks(CurIDX) = tracks_img_to_world(FinalTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
end

EmptyTracksHeader = TracksHeader;

TracksHeader.Tracks = FinalTracks;
TracksHeader.count = length(TracksHeader.Tracks);
F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
save_mif_tracks(TracksHeader, [F '.tck']);
clear TracksHeader;
[~,~,~] = mkdir(fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois']));
for I = 1:length(CurLabels.shortlabels)
	T = EmptyTracksHeader;
	T.Tracks = FinalTracks((FinalStartRegions == I) | (FinalEndRegions == I));
	F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois'], [CurLabels.shortlabels{I} '.tck']);
	save_mif_tracks(T, F);
	clear T F;
end

TracksHeader = EmptyTracksHeader;
TracksHeader.Tracks = FinalTracks;

F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
save_trackvis_tracks_nii([F '.trk'], TracksHeader, WMMaskNII);
clear T2;

[~,~,~] = mkdir(fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois']));
for I = 1:length(CurLabels.shortlabels)
	T = EmptyTracksHeader;
	T.Tracks = FinalTracks((FinalStartRegions == I) | (FinalEndRegions == I));
	F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois'], [CurLabels.shortlabels{I} '.trk']);
	save_trackvis_tracks_nii(F, T, WMMaskNII);
	clear T F;
end

SeedSizeMatrix = bsxfun(@plus, SeedSizes(:), SeedSizes(:)');
SeedSizeMatrix = SeedSizeMatrix + double(SeedSizeMatrix == 0);
SizeWeightedA = CountA ./ SeedSizeMatrix;
WeightedA = WeightedA .* 2 ./ SeedSizeMatrix;
T = CountA;
T(T == 0) = 1;
LengthA = LengthA ./ T;
clear T;

save(fullfile(Subject, ['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalStartRegions', 'FinalEndRegions', 'WeightedA', 'CountA', 'SizeWeightedA', 'LengthA', 'SeedSizes');

AcceptedSeedsNII.img = int32(permute(flipdim(AcceptedSeedsIMG, 1), [2 1 3]));
save_nii(AcceptedSeedsNII, fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_seeds_' num2str(Curvature) '.nii.gz']));
%end


function [ID, JD, KD] = first_last_coords_to_idx(S)

ID = S([1 end], 2);
JD = S([1 end], 1);
KD = S([1 end], 3);

ID = ID(:)';
JD = JD(:)';
KD = KD(:)';

function [T] = first_last_coords(S)

T = S([1 end], :);

function [StartPoint, EndPoint] = streamline_start_end_points(S)

StartPoint = S(1, :);
EndPoint = S(end, :);
%T = S([1 end], :);


function [StartVector, EndVector] = streamline_start_end_vectors(S)

if(size(S) > 2)
	StartVector = S(1, :) - S(2, :);
	EndVector = S(end, :) - S(end - 1, :);
else
	StartVector = [0, 0, 0];
	EndVector = [0, 0, 0];
end


function [M] = tracks_anterior_cc_posterior_trajectory(S, SeedIMG)

%M = false(size(S));

% 253 is the central CC label
I = find(SeedIMG == 253);
[ID, ~, ~] = ind2sub(size(SeedIMG), I);
CentralCCCentroid = mean(ID);

SZSlice = int32(size(SeedIMG, 1) * size(SeedIMG, 2));
SZRows = int32(size(SeedIMG, 1));

T = cat(1, S{:});
TracksSZ = cellfun('size', S, 1);

PastCCCentroid = (T(:, 2) > CentralCCCentroid);
PastCCCentroid = mat2cell_vec(PastCCCentroid, TracksSZ);
PastCCCentroid = cellfun(@any, PastCCCentroid);

T = int32(round(T));
IDX = SZSlice * (T(:, 3) - 1) + SZRows * (T(:, 2) - 1) + T(:, 1);

SeedValues = SeedIMG(IDX);
SeedValuesCC = (SeedValues == 255);
clear SeedValues;
SeedValuesCC = mat2cell_vec(SeedValuesCC, TracksSZ);
SeedValuesCC = cellfun(@any, SeedValuesCC);

M = SeedValuesCC & PastCCCentroid;

function [FornixTracks, AnyInterhemispheric, AnyInterhemisphericNoCC, AnyInterhemisphericAnteriorCC, AnyInterhemisphericPosteriorCC] = tracks_through_fornix(S, SeedIMG)

% the labels for the left cerebral white matter and right cerebral white
% matter are 2 and 41, respectively, so if we find a [2 41] or a [41 2] the
% track is marked as 1

FornixTracks = false(size(S));
%InterhemisphericNoCC = false(size(S));
AnyInterhemispheric = false(size(S));
AnyInterhemisphericNoCC = false(size(S));
AnyInterhemisphericAnteriorCC = false(size(S));
AnyInterhemisphericPosteriorCC = false(size(S));
T = cat(1, S{:});
TracksSZ = cellfun(@numrows, S);

SeedValues = interp3(double(SeedIMG), double(T(:, 1)), double(T(:, 2)), double(T(:, 3)), 'nearest');
SeedValues = mat2cell_vec(SeedValues, TracksSZ);
clear T;
CCLabels = [251 252 253 254 255];
for z = 1:length(S)
	LWM = (SeedValues{z} == 2);
	RWM = (SeedValues{z} == 41);

	if(any(LWM) && any(RWM))
		AnyInterhemispheric(z) = 1;
	end

	if(any(SeedValues{z} == 251))
		AnyInterhemisphericPosteriorCC(z) = 1;
	end
	
	if(any(SeedValues{z} == 255))
		AnyInterhemisphericAnteriorCC(z) = 1;
	end
	
	if(AnyInterhemispheric(z) == 1)
		if(all(~ismember(SeedValues{z}, CCLabels)))
			AnyInterhemisphericNoCC(z) = 1;
		end
	end
	
	
	if(AnyInterhemispheric(z) == 1)
		if(all(~ismember(SeedValues{z}, CCLabels)))
			AnyInterhemisphericNoCC(z) = 1;
		end
	end
	if(any(...
			(RWM(2:end) & LWM(1:end - 1)) | ...
			(LWM(2:end) & RWM(1:end - 1)) ...
			))
		%keyboard;
		FornixTracks(z) = 1;
	end
end

function M = tracks_no_wm_mask(S, SeedIMG)

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

%NotWM = SeedValues < 0.5;%~((SeedValues == 2) | (SeedValues == 41));%ismember(SeedValues, [2 41]);

%keyboard;
%clear SeedValues;


%M2 = cellfun(@all, SeedValues);
% for k = 1:length(NotWM)
% 	if(all(NotWM{k}))
% 		M(CurIDX(k)) = 1;
% 	end
% end
% 
% end
% T = cat(1, S{:});
% TracksSZ = cellfun('size', S, 1);
% 
% SeedValues = interp3_linear_fast(double(SeedIMG), double(T(:, 1)), double(T(:, 2)), double(T(:, 3)), 'nearest');
% %SeedValues = mat2cell_vec(SeedValues, TracksSZ);
% NotWM = ~((SeedValues == 2) | (SeedValues == 41));%ismember(SeedValues, [2 41]);
% NotWM = mat2cell_vec(NotWM, TracksSZ);
% clear SeedValues;
% clear T;
% 
% for z = 1:length(S)
% 	if(all(NotWM{z}))
% 		M(z) = 1;
% 	end
% end
