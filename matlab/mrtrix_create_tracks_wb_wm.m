function [varargout] = mrtrix_create_tracks_wb_wm(Subject, MrtrixMethod, SeedType, Curvature)

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
% PARAMETERS
%	Subject (string): the ID of the subject, must be a directory
%	MrtrixMethod (string): tractography method 'dt_stream', 'sd_stream', or 'sd_prob'
%	SeedType (string): the parcellation scheme to use
%		Freesurfer defaults are 'aparc', 'aparc.a2005s', 'aparc.a2009s'
%		For others see the available files
%	Curvature [1]: -curvature radius parameter for streamtrack
%		Hagmann's "Mapping Structural Core of the Cerebral Cortex" uses 30 degrees per 1mm, which is 1.9
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

%MrtrixSubjDir = mrtrix_create_tracks_mrtrixdir;

MrtrixSubjDir = 'mrtrix2DWISpaceBReorient';

[CurLabels, IncludeFile, ExcludeFile] = mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature);

WMMaskNII = load_nii(fullfile(MrtrixSubjDir, Subject, 'freesurfer_wm'));
WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);

WMMaskIMG = (WMMaskIMG > 0);

% this creates too many tracks for high res volumes
% create two tracks per voxel of 0.75 x 0.75 x 0.75
WMMaskNIIVoxVolume = prod(WMMaskNII.hdr.dime.pixdim(2:4));

DesiredVolume = 0.7 * 0.7 * 0.7;
%DesiredVolume = 1;
TotalNumTracks = round(WMMaskNIIVoxVolume / DesiredVolume * sum(WMMaskIMG(:)));
clear DesiredVolume;
%TotalNumTracks = 10000;
disp(['Total tracks: ' num2str(TotalNumTracks)]);
NumTracksSoFar = 0;
FinalTracks = cell(1, TotalNumTracks);
FinalArcLengths = zeros(1, TotalNumTracks);
FinalStartRegions = zeros(1, TotalNumTracks);
FinalEndRegions = zeros(1, TotalNumTracks);
TracksHeader = [];

OriginalSeedNII = load_nii(fullfile(MrtrixSubjDir, '..', 'RegisterFreesurferToMrtrix', Subject, CurLabels.SeedFile));
OriginalSeedIMG = uint16(flipdim(permute(OriginalSeedNII.img, [2 1 3]), 1));
OriginalSeedNII.img = [];

SeedNII = load_nii(IncludeFile);
SeedIMG = uint16(flipdim(permute(SeedNII.img, [2 1 3]), 1));
SeedNII.img = [];

FANII = load_nii(fullfile(MrtrixSubjDir, Subject, 'fa'));
FAIMG = flipdim(permute(FANII.img, [2 1 3]), 1);
FANII.img = [];

T = imdilate(WMMaskIMG, ones(3, 3, 3));
SeedBoundaries = SeedIMG(T);
SeedBoundaries = SeedBoundaries(SeedBoundaries > 0);
clear T;

% disp('Saving seeds');
% SeedOutputDir = fullfile(MrtrixSubjDir, Subject, ['seeds_' SeedType]);
% [~,~,~] = mkdir(SeedOutputDir);
% for z = 1:length(CurLabels.values)
%  	CurSeedIMG = uint8(SeedIMG == z);
%  	CurSeedNII = SeedNII;
%  	CurSeedNII.img = permute(flipdim(CurSeedIMG, 1), [2 1 3]);
%  	save_nii(CurSeedNII, fullfile(MrtrixSubjDir, Subject, ['seeds_' SeedType], [CurLabels.shortlabels{z} '.nii.gz']));
%  	clear CurSeedIMG CurSeedNII;
% end
% disp('Finished saving seeds');

NumRegions = length(CurLabels.values);

SeedVolumes = histc(SeedIMG(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;
SeedSurfaceAreas = histc(SeedBoundaries(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;
clear SeedBoundaries;

AcceptedSeedsNII = WMMaskNII;
AcceptedSeedsIMG = zeros(size(SeedIMG), 'int32');

WeightedA = zeros(length(CurLabels.values));
LengthA = zeros(length(CurLabels.values));
CountA = zeros(length(CurLabels.values));

MrtrixBinDir = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9-accepted-only-cout-avx2')];
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
		fullfile(MrtrixBinDir, 'lib'), ...
		fullfile(MrtrixBinDir, 'bin', 'streamtrack'), ...
		fullfile(MrtrixSubjDir, Subject, 'dwi.nii.gz'), ...
		TrackFile, ...
		fullfile(MrtrixSubjDir, Subject, 'freesurfer_wm.nii.gz'), ...
		ExcludeFile, ...
		IncludeFile, ...
		fullfile(MrtrixSubjDir, Subject, 'mask.nii.gz'), NumTracks, NumTracks * 100, ...
		MrtrixStep, ...
		fullfile(MrtrixSubjDir, Subject, 'grad.b'), ...
		Curvature, ...
		10, ...
		LogFile);
		case {'sd_stream', 'sd_prob', 'sd_runge2', 'sd_runge4'}
	%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
	CommandString = sprintf('LD_LIBRARY_PATH=%s; %s %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -stop -curvature %f -minlength %f > %s;', ...
		fullfile(MrtrixBinDir, 'lib'), ...
		fullfile(MrtrixBinDir, 'bin', 'streamtrack'), ...
		upper(MrtrixMethod), ...
		fullfile(MrtrixSubjDir, Subject, 'CSD.nii.gz'), ...
		TrackFile, ...
		fullfile(MrtrixSubjDir, Subject, 'freesurfer_wm.nii.gz'), ...
		ExcludeFile, ...
		IncludeFile, ...
		fullfile(MrtrixSubjDir, Subject, 'mask.nii.gz'), NumTracks, NumTracks * 100, ...
		MrtrixStep, ...
		fullfile(MrtrixSubjDir, Subject, 'grad.b'), ...
		Curvature, ...
		10, ...
		LogFile);
	end
	disp(CommandString);
	disp('Generating tracks');
	system(CommandString);
	clear CommandString;
	disp('Finished generating tracks');
	MIFTracks = load_mif_tracks(TrackFile);
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
	clear LogFile;
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
		clear T M;
		clear RightIDX CurIDX;
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
	clear RightIDX CurIDX;
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

	clear tokens Codes;
	
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
	clear UniqueCodes TracksToAddMask TracksToAddIDX TracksOutsideIMG TracksSZ StartRegion EndRegion NumBefore* IDX CurTracksSZ AcceptedTrackMask ArcLengths;
	%keyboard;
end

% make sure that the start regions are less than the end regions
T = sort([FinalStartRegions(:), FinalEndRegions(:)], 2, 'ascend');
[FinalStartRegions, FinalEndRegions] = deal(T(:, 1)', T(:, 2)');
clear T;

[T, I] = sortrows([FinalStartRegions(:), FinalEndRegions(:)], [1, 2]);
[FinalStartRegions, FinalEndRegions] = deal(T(:, 1)', T(:, 2)');

FinalTracks = FinalTracks(I);
FinalArcLengths = FinalArcLengths(I);
% reorder tracks so that they are in sorted order of edges
%sort([FinalStartRegions(:), FinalEndRegions(:)], 2, 'ascend');

% TRACKS ARE IN IMAGE SPACE
%save(fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalTracks', 'FinalArcLengths', 'FinalStartRegions', 'FinalEndRegions', '-v7.3');

% no longer are we saving track files automatically, do it on request
NumFinalTracks = numel(FinalTracks);
FinalTracksSZ = cellfun('size', FinalTracks, 1);
% block to save RAM
%FinalTracks = tracks_img_to_world(FinalTracks, WMMaskNII, size(WMMaskIMG));

FinalTracksWorld = cell(size(FinalTracks));
for z = 1:BlockSize:NumFinalTracks
	RightIDX = min(z + BlockSize - 1, NumFinalTracks);
 	CurIDX = (z:RightIDX);
 
 	FinalTracksWorld(CurIDX) = tracks_img_to_world(FinalTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
end
clear RightIDX CurIDX;

FATracks = tracks_world_to_img(FinalTracksWorld, FANII, size(FAIMG));
FATracks = cell2mat(FATracks);

FAOfTracks = interp3(double(FAIMG), FATracks(:, 1), FATracks(:, 2), FATracks(:, 3));
clear FATracks;

FAOfTracks = mat2cell_vec(FAOfTracks, FinalTracksSZ(:));

MeanFAOfTracks = cellfun(@mean, FAOfTracks);

MeanFAA = zeros(size(CountA));

%RegionIDX = sub2ind(size(CountA), FinalStartRegions, FinalEndRegions);
%

for RegionOne = 1:NumRegions - 1
	for RegionTwo = RegionOne + 1:NumRegions
		I = (FinalStartRegions == RegionOne) & (FinalEndRegions == RegionTwo);
		if(any(I))
			MeanFAA(RegionOne, RegionTwo) = mean(MeanFAOfTracks(I));
		end
	end
end
MeanFAA = MeanFAA + MeanFAA';
%keyboard;
% 
% %%
% 
% SR = 1;
% SC = 2;
% 
% subplot(SR, SC, 1);
% cla;
% s1 = slice(double(SeedIMG), floor(size(SeedIMG, 2) / 2) + 0, floor(size(SeedIMG, 1) / 2) + 0, floor(size(SeedIMG, 3) / 2) + 0);
% set(s1, 'EdgeColor', 'none');
% colormap gray;
% axis ij equal;
% H1 = streamline(FinalTracks(1:5:end));
% subplot(SR, SC, 2);
% cla;
% s2 = slice(double(FAIMG), floor(size(FAIMG, 2) / 2) + 0, floor(size(FAIMG, 1) / 2) + 0, floor(size(FAIMG, 3) / 2) + 0);
% set(s2, 'EdgeColor', 'none');
% axis ij equal;
% H2 = streamline(FATracks(1:5:end));

%%
%
% EmptyTracksHeader = TracksHeader;
% 
% TracksHeader.Tracks = FinalTracksWorld;
% TracksHeader.count = length(TracksHeader.Tracks);
% F = fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
% save_mif_tracks(TracksHeader, [F '.tck']);
% clear TracksHeader;
% [~,~,~] = mkdir(fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois']));
% for I = 1:length(CurLabels.shortlabels)
% 	T = EmptyTracksHeader;
% 	T.Tracks = FinalTracksWorld((FinalStartRegions == I) | (FinalEndRegions == I));
% 	F = fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois'], [CurLabels.shortlabels{I} '.tck']);
% 	save_mif_tracks(T, F);
% 	clear T F;
% end
% 
% TracksHeader = EmptyTracksHeader;
% TracksHeader.Tracks = FinalTracksWorld;
% 
% F = fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
% save_trackvis_tracks_nii([F '.trk'], TracksHeader, WMMaskNII);
% clear T2;
% 
% [~,~,~] = mkdir(fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois']));
% for I = 1:length(CurLabels.shortlabels)
% 	T = EmptyTracksHeader;
% 	T.Tracks = FinalTracks((FinalStartRegions == I) | (FinalEndRegions == I));
% 	F = fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois'], [CurLabels.shortlabels{I} '.trk']);
% 	save_trackvis_tracks_nii(F, T, WMMaskNII);
% 	clear T F;
% end

SeedVolumeSumMatrix = bsxfun(@plus, SeedVolumes(:), SeedVolumes(:)');
SeedVolumeSumMatrix = SeedVolumeSumMatrix + double(SeedVolumeSumMatrix == 0);
SeedVolumeProdMatrix = bsxfun(@times, SeedVolumes(:), SeedVolumes(:)');
SeedVolumeProdMatrix = SeedVolumeProdMatrix + double(SeedVolumeProdMatrix == 0);

SeedSurfaceAreaSumMatrix = bsxfun(@plus, SeedSurfaceAreas(:), SeedSurfaceAreas(:)');
SeedSurfaceAreaSumMatrix = SeedSurfaceAreaSumMatrix + double(SeedSurfaceAreaSumMatrix == 0);
SeedSurfaceAreaProdMatrix = bsxfun(@times, SeedSurfaceAreas(:), SeedSurfaceAreas(:)');
SeedSurfaceAreaProdMatrix = SeedSurfaceAreaProdMatrix + double(SeedSurfaceAreaProdMatrix == 0);

VolumeSumCountA = CountA ./ SeedVolumeSumMatrix;
VolumeProdCountA = CountA ./ SeedVolumeProdMatrix;

SurfaceAreaSumCountA = CountA ./ SeedSurfaceAreaSumMatrix;
SurfaceAreaProdCountA = CountA ./ SeedSurfaceAreaProdMatrix;

VolumeSumWeightedA = WeightedA .* 2 ./ SeedVolumeSumMatrix;
VolumeProdWeightedA = WeightedA .* 2 ./ SeedVolumeProdMatrix;

SurfaceAreaSumWeightedA = WeightedA .* 2 ./ SeedSurfaceAreaSumMatrix;
SurfaceAreaProdWeightedA = WeightedA .* 2 ./ SeedSurfaceAreaProdMatrix;

T = CountA;
T(T == 0) = 1;
LengthA = LengthA ./ T;
clear T;

save(fullfile(MrtrixSubjDir, Subject, ['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), ...
	'FinalStartRegions', 'FinalEndRegions', ...
	'SeedVolumes', 'SeedSurfaceAreas', ...
	'CountA', 'LengthA', 'MeanFAA', ...
	'VolumeSumCountA', 'VolumeProdCountA', ...
	'SurfaceAreaSumCountA', 'SurfaceAreaProdCountA', ...
	'VolumeSumWeightedA', 'VolumeProdWeightedA', ...
	'SurfaceAreaSumWeightedA', 'SurfaceAreaProdWeightedA', ...
	'-v7.3');

AcceptedSeedsNII.img = int32(permute(flipdim(AcceptedSeedsIMG, 1), [2 1 3]));
save_nii(AcceptedSeedsNII, fullfile(MrtrixSubjDir, Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_seeds_' num2str(Curvature) '.nii.gz']));
