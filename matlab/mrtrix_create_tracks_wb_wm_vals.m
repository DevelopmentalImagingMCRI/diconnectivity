function [varargout] = mrtrix_create_tracks_wb_wm_vals(Subject, MrtrixMethod, SeedType, Curvature)

% mrtrix_create_tracks_wb_wm(Subject, MrtrixMethod, SeedType, Curvature)
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

[CurLabels, IncludeFile, ExcludeFile] = mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature);

WMMaskNII = load_nii(fullfile(Subject, 'freesurfer_wm'));
WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);

% this creates too many tracks for high res volumes
% create two tracks per voxel of 0.75 x 0.75 x 0.75
WMMaskNIIVoxVolume = prod(WMMaskNII.hdr.dime.pixdim(2:4));

DesiredVolume = 0.7 * 0.7 * 0.7;
%DesiredVolume = 1;
TotalNumTracks = round(WMMaskNIIVoxVolume / DesiredVolume * sum(WMMaskIMG(:) > 0));
clear DesiredVolume;
TotalNumTracks = 10;
disp(['Total tracks: ' num2str(TotalNumTracks)]);
NumTracksSoFar = 0;
FinalTracks = cell(1, TotalNumTracks);
FinalArcLengths = zeros(1, TotalNumTracks);
FinalStartRegions = zeros(1, TotalNumTracks);
FinalEndRegions = zeros(1, TotalNumTracks);
TracksHeader = [];

OriginalSeedNII = load_nii(fullfile('..', 'RegisterFreesurferToMrtrix', Subject, CurLabels.SeedFile));
OriginalSeedIMG = flipdim(permute(OriginalSeedNII.img, [2 1 3]), 1);
SeedNII = load_nii(IncludeFile);
SeedIMG = flipdim(permute(SeedNII.img, [2 1 3]), 1);

% disp('Saving seeds');
% SeedOutputDir = fullfile(Subject, ['seeds_' SeedType]);
% [~,~,~] = mkdir(SeedOutputDir);
% for z = 1:length(CurLabels.values)
%  	CurSeedIMG = uint8(SeedIMG == z);
%  	CurSeedNII = SeedNII;
%  	CurSeedNII.img = permute(flipdim(CurSeedIMG, 1), [2 1 3]);
%  	save_nii(CurSeedNII, fullfile(Subject, ['seeds_' SeedType], [CurLabels.shortlabels{z} '.nii.gz']));
%  	clear CurSeedIMG CurSeedNII;
% end
% disp('Finished saving seeds');

SeedSizes = histc(SeedIMG(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;

AcceptedSeedsNII = WMMaskNII;
AcceptedSeedsIMG = zeros(size(SeedIMG), 'int32');

WeightedA = zeros(length(CurLabels.values));
LengthA = zeros(length(CurLabels.values));
CountA = zeros(length(CurLabels.values));

MrtrixDirectory = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9-accepted-only-cout-prob')];
BlockSize = 50000;

while(NumTracksSoFar < TotalNumTracks)
	MrtrixStep = 1;
	NumTracks = max(min((TotalNumTracks - NumTracksSoFar) * 2, 2000000), 1000);
	TrackFile = [tempname '.tck'];
	LogFile = [tempname '.log'];
	switch(lower(MrtrixMethod))
		case 'dt_stream'
	%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
	CommandString = sprintf('LD_LIBRARY_PATH=%s; %s DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -curvature %f -stop -minlength %f > %s;', ...
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
		10, ...
		LogFile);
		case {'sd_stream', 'sd_prob', 'sd_runge2', 'sd_runge4'}
	%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
	CommandString = sprintf('LD_LIBRARY_PATH=%s; %s %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -stop -curvature %f -minlength %f > %s 2>streamtrack.out;', ...
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
		10, ...
		LogFile);
	end
	disp(CommandString);
	disp('Generating tracks');
	system(CommandString);
	disp('Finished generating tracks');
	MIFTracks = load_mifvals_tracks(TrackFile);
	
	for z = 1:length(MIFTracks.Vals)
		MIFTracks.Vals{z}(end) = MIFTracks.Vals{z}(end - 1);
		T = [MIFTracks.Vals{z}(1), MIFTracks.Vals{z}(:)', MIFTracks.Vals{z}(end)];
		MIFTracks.Vals{z} = conv(T, [0.25 0.5 0.25], 'valid');
		clear T;
	end
	keyboard;
	%%
	clf;
	%M = cellfun(@max, MIFTracks.Vals);
	CMAPX = linspace(0, 1, 256);
	CMAP = jet(256);
	
	for z = 1:1000%length(MIFTracks.Tracks)
		MIFTracks.Vals{z} = min(MIFTracks.Vals{z}, 1);
		MIFTracks.Vals{z} = max(MIFTracks.Vals{z}, 0);
		M = mean(MIFTracks.Vals{z});
		C = zeros(1, 3);
		C(1) = interp1q(CMAPX, CMAP(:, 1), M);
		C(2) = interp1q(CMAPX, CMAP(:, 2), M);
		C(3) = interp1q(CMAPX, CMAP(:, 3), M);
		L = streamline(MIFTracks.Tracks(z));
		set(L, 'Color', C);
% 		for k = 1:size(MIFTracks.Tracks{z}, 1) - 1
% 			C = zeros(1, 3);
% 			C(1) = interp1q(CMAPX, CMAP(:, 1), MIFTracks.Vals{z}(k));
% 			C(2) = interp1q(CMAPX, CMAP(:, 2), MIFTracks.Vals{z}(k));
% 			C(3) = interp1q(CMAPX, CMAP(:, 3), MIFTracks.Vals{z}(k));
% 			C(isnan(C)) = 0;
% 			L = line([MIFTracks.Tracks{z}(k, 1) MIFTracks.Tracks{z}(k + 1, 1)], ... 
% 				[MIFTracks.Tracks{z}(k, 2) MIFTracks.Tracks{z}(k + 1, 2)], ...
% 				[MIFTracks.Tracks{z}(k, 3) MIFTracks.Tracks{z}(k + 1, 3)]);
% 			set(L, 'Color', C);
% 		end
	end
	set(gca, 'Color', 'k');
	axis equal;
	%%
	%copyfile(TrackFile, 'tmp.tck');
	delete(TrackFile);
	CurTracks = MIFTracks.Tracks;
	if(isempty(TracksHeader))
		TracksHeader = MIFTracks;
		TracksHeader.Tracks = [];
	end
	MIFTracks.Tracks = [];

	FID = fopen(LogFile, 'r');
	tokens = textscan(FID, '%s %f %f %f %d %d %d %d %d');
	fclose(FID);
	delete(LogFile);
	for z = 1:length(tokens)
		tokens{z} = tokens{z}(1:MIFTracks.count);
	end

	clear MIFTracks;

	TracksSZ = cellfun('size', CurTracks, 1);
	CurNumTracks = numel(CurTracks);
	ArcLengths = TracksSZ * MrtrixStep;

	% convert the tracks from world coordinates (mrtrix) to voxel
	% coordinates
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = (z:RightIDX);

		CurTracks(CurIDX) = tracks_world_to_img(CurTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
	end

	AcceptedTrackMask = strcmp(tokens{1}, 'accept:');
	WMMaskIMGSZ = size(WMMaskIMG);
	TracksOutsideIMG = false(size(AcceptedTrackMask));
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = (z:RightIDX);
		CurTracksSZ = TracksSZ(CurIDX);
		
		T = cat(1, CurTracks{CurIDX});
		M = any(T < 1, 2) | ...
			T(:, 1) > WMMaskIMGSZ(2) | ...
			T(:, 2) > WMMaskIMGSZ(1) | ...
			T(:, 3) > WMMaskIMGSZ(3);
		M = mat2cell_vec(M, int32(CurTracksSZ));
		TracksOutsideIMG(CurIDX) = cellfun(@any, M);
		clear T;
	end
	
	AcceptedTrackMask = AcceptedTrackMask & ~TracksOutsideIMG;
	disp([num2str(sum(AcceptedTrackMask)) ' tracks left after insideimg']);
	
	
	disp([num2str(sum(AcceptedTrackMask)) ' mrtrix accepted tracks']);
	disp([num2str(sum(~AcceptedTrackMask)) ' mrtrix rejected tracks']);

	StartRegion = tokens{8};
	EndRegion = tokens{9};

	AcceptedTrackMask(StartRegion == EndRegion) = 0;
	disp([num2str(sum(AcceptedTrackMask)) ' tracks left after same region cull']);
	NumBeforeNoWM = sum(AcceptedTrackMask);

	CurNumTracks = sum(AcceptedTrackMask);
	IDX = find(AcceptedTrackMask);
	RejectedTracksNoWM = false(CurNumTracks, 1);
	for z = 1:BlockSize:CurNumTracks
		RightIDX = min(z + BlockSize - 1, CurNumTracks);
		CurIDX = z:RightIDX;
		RejectedTracksNoWM(CurIDX) = mrtrix_create_tracks_no_wm_mask(CurTracks(IDX(CurIDX)), WMMaskIMG);
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
		RejectedAnteriorCCPosterior(CurIDX) = mrtrix_create_tracks_anterior_cc_posterior_trajectory(CurTracks(IDX(CurIDX)), OriginalSeedIMG);
	end
	AcceptedTrackMask(AcceptedTrackMask) = ~RejectedAnteriorCCPosterior;

	clear RejectedAnteriorCCPosterior;
	disp([num2str(NumBeforeNoWM - NumBeforeAnteriorPosterior) ' culled by no wm']);
	disp([num2str(NumBeforeAnteriorPosterior - sum(AcceptedTrackMask)) ' culled by anterior posterior']);
	%$#AcceptedTracks = CurTracks(AcceptedTrackMask);

	Seeds = cat(2, tokens{2:4});

	NumTracksToAdd = min(sum(AcceptedTrackMask), TotalNumTracks - NumTracksSoFar);
	TracksToAddIDX = find(AcceptedTrackMask, NumTracksToAdd, 'first');
	TracksToAddMask = false(size(AcceptedTrackMask));
	TracksToAddMask(TracksToAddIDX) = 1;
	TracksToAdd = CurTracks(TracksToAddMask);

	ToAddSeeds = Seeds(TracksToAddIDX, :);

	clear Seeds;
	Codes = tokens{5};

	UniqueCodes = unique(Codes);
	for z = 1:length(UniqueCodes)
		disp([num2str(sum(Codes == UniqueCodes(z))) ' got code ' num2str(UniqueCodes(z))]);
	end

	clear tokens;

	disp([num2str(sum(~AcceptedTrackMask)) ' rejected tracks']);
	disp([num2str(sum(AcceptedTrackMask)) ' accepted tracks']);
	disp([num2str(sum(~TracksToAddMask)) ' rejected tracks along with those that arent added']);

	ToAddArcLengths = ArcLengths(TracksToAddMask);
	ToAddStartRegions = StartRegion(TracksToAddMask);
	ToAddEndRegions = EndRegion(TracksToAddMask);
	clear RejectedTracks;

	ToAddSeedsVox = tracks_world_to_img({ToAddSeeds}, WMMaskNII, size(WMMaskIMG));
	ToAddSeedsVox = ToAddSeedsVox{1};

	I = round(ToAddSeedsVox(:, 2)) + size(AcceptedSeedsIMG, 1) * round(ToAddSeedsVox(:, 1)) + size(AcceptedSeedsIMG, 1) * size(AcceptedSeedsIMG, 2) * round(ToAddSeedsVox(:, 3));
	AddToIMG = histc(I, 1:numel(AcceptedSeedsIMG));
	AcceptedSeedsIMG = AcceptedSeedsIMG + int32(reshape(AddToIMG, size(AcceptedSeedsIMG)));
	clear I AddToIMG AcceptedSeedsVox;

	% fast way of adding the counts to the connectivity matrices
	% get the linear indices of the regions in the connectivity
	% matrices

	CurRegionsI = sub2ind(size(CountA), ToAddStartRegions, ToAddEndRegions);

	[CurRegionsISorted, CurRegionsISortedIIDX] = sort(CurRegionsI);

	% get the sizes and indices of 
	FirstIndices = [1; find(diff(CurRegionsISorted) > 0) + 1];
	CurSZ = diff([FirstIndices; length(CurRegionsI) + 1]);

	LengthAValues = mat2cell_vec(ToAddArcLengths(CurRegionsISortedIIDX), CurSZ);
	WeightedAValues = mat2cell_vec(1 ./ ToAddArcLengths(CurRegionsISortedIIDX), CurSZ);

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

	FinalTracks(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd) = TracksToAdd;
	FinalArcLengths(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd) = ToAddArcLengths;

	FinalStartRegions(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd) = ToAddStartRegions;
	FinalEndRegions(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd) = ToAddEndRegions;
	NumTracksSoFar = NumTracksSoFar + NumTracksToAdd;

	clear CurTracks SeedValues RejectedTracks TracksToAdd NumTracksToAdd;
	clear ToAdd*;
	disp(['Number of tracks generated so far: ' num2str(NumTracksSoFar) ' of ' num2str(TotalNumTracks)]);
end

save(fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalTracks', 'FinalArcLengths', 'FinalStartRegions', 'FinalEndRegions', '-v7.3');

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

save(fullfile(Subject, ['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalStartRegions', 'FinalEndRegions', 'WeightedA', 'CountA', 'SizeWeightedA', 'LengthA', 'SeedSizes', '-v7.3');

AcceptedSeedsNII.img = int32(permute(flipdim(AcceptedSeedsIMG, 1), [2 1 3]));
save_nii(AcceptedSeedsNII, fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_seeds_' num2str(Curvature) '.nii.gz']));
