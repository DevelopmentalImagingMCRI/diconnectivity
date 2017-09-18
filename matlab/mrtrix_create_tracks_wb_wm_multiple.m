function [varargout] = mrtrix_create_tracks_wb_wm_multiple(Subject, MrtrixMethod, SeedType, Curvature, NumIters)

% mrtrix_create_tracks_wb_wm_multiple(Subject, MrtrixMethod, SeedType, Curvature)
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
TotalNumTracks = round(WMMaskNIIVoxVolume / DesiredVolume * sum(WMMaskIMG(:) > 0));
clear DesiredVolume;
%TotalNumTracks = 5000;

OriginalSeedNII = load_nii(fullfile('..', 'RegisterFreesurferToMrtrix', Subject, CurLabels.SeedFile));
OriginalSeedIMG = flipdim(permute(OriginalSeedNII.img, [2 1 3]), 1);
SeedNII = load_nii(IncludeFile);
SeedIMG = flipdim(permute(SeedNII.img, [2 1 3]), 1);

disp('Saving seeds');
SeedOutputDir = fullfile(Subject, ['seeds_' SeedType]);
[~,~,~] = mkdir(SeedOutputDir);
for z = 1:length(CurLabels.values)
 	CurSeedIMG = uint8(SeedIMG == z);
 	CurSeedNII = SeedNII;
 	CurSeedNII.img = permute(flipdim(CurSeedIMG, 1), [2 1 3]);
 	save_nii(CurSeedNII, fullfile(Subject, ['seeds_' SeedType], [CurLabels.shortlabels{z} '.nii.gz']));
 	clear CurSeedIMG CurSeedNII;
end
disp('Finished saving seeds');

SeedSizes = histc(SeedIMG(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;

%AcceptedSeedsNII = WMMaskNII;
%AcceptedSeedsIMG = zeros(size(SeedIMG), 'int32');

NumROIs = length(CurLabels.values);
WeightedA = zeros([NumROIs, NumROIs, NumIters]);
SizeWeightedA = zeros([NumROIs, NumROIs, NumIters]);
LengthA = zeros([NumROIs, NumROIs, NumIters]);
CountA = zeros([NumROIs, NumROIs, NumIters]);

MrtrixDirectory = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9-accepted-only-cout')];
BlockSize = 50000;

%FinalTracks = cell(TotalNumTracks, NumIters);
%FinalArcLengths = zeros(TotalNumTracks, NumIters);
FinalStartRegions = zeros(TotalNumTracks, NumIters);
FinalEndRegions = zeros(TotalNumTracks, NumIters);
FinalArcLengths = zeros(TotalNumTracks, NumIters);
TracksHeader = [];

for CurIter = 1:NumIters
	disp(['Iteration ' num2str(CurIter)]);
	NumTracksSoFar = 0;
	while(NumTracksSoFar < TotalNumTracks)
		MrtrixStep = 1;
		NumTracks = max(min((TotalNumTracks - NumTracksSoFar) * 2, 2000000), 1000);
		TrackFile = [tempname '.tck'];
		ResultFile = tempname;
		switch(lower(MrtrixMethod))
			case 'dt_stream'
		%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
		CommandString = sprintf('LD_LIBRARY_PATH=%s; %s DT_STREAM %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -curvature %f -stop -minlength %f 2>/dev/null > %s;', ...
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
			ResultFile);
			case {'sd_stream', 'sd_prob', 'sd_runge2', 'sd_runge4'}
		%CommandString = sprintf('LD_LIBRARY_PATH=/usr/local/mrtrix-0.2.9/lib; streamtrack %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step 1 -grad %s -curvature 4 -stop 2>/dev/null;', ...
		CommandString = sprintf('LD_LIBRARY_PATH=%s; %s %s %s %s -quiet -seed %s -exclude %s -include %s -mask %s -number %d -maxnum %d -step %f -grad %s -stop -curvature %f -minlength %f > %s;', ...
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
			ResultFile);
		end
		disp(CommandString);
		disp('Generating tracks');
		system(CommandString);
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
		FID = fopen(ResultFile, 'r');
		tokens = textscan(FID, '%s %f %f %f %d %d %d %d %d');
		fclose(FID);
		delete(ResultFile);
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
				
		NumTracksToAdd = min(sum(AcceptedTrackMask), TotalNumTracks - NumTracksSoFar);
		ToAddTracksIDX = find(AcceptedTrackMask, NumTracksToAdd, 'first');
		ToAddTrackMask = false(size(AcceptedTrackMask));
		ToAddTrackMask(ToAddTracksIDX) = 1;
		
		%TracksToAdd = CurTracks(ToAddTrackMask);
		
		%Seeds = cat(2, tokens{2:4});
		%ToAddSeeds = Seeds(ToAddTracksIDX, :);
				
		clear Seeds;
		Codes = tokens{5};

		UniqueCodes = unique(Codes);
		for z = 1:length(UniqueCodes)
			disp([num2str(sum(Codes == UniqueCodes(z))) ' got code ' num2str(UniqueCodes(z))]);
		end
		clear tokens;
		
		disp([num2str(sum(~AcceptedTrackMask)) ' rejected tracks']);
		disp([num2str(sum(AcceptedTrackMask)) ' accepted tracks']);
		disp([num2str(sum(~ToAddTrackMask)) ' rejected tracks along with those that arent added']);
		
		ToAddArcLengths = ArcLengths(ToAddTrackMask);
		ToAddStartRegions = StartRegion(ToAddTrackMask);
		ToAddEndRegions = EndRegion(ToAddTrackMask);
		
		clear RejectedTracks;
		
		% fast way of adding the counts to the connectivity matrices
		% get the linear indices of the regions in the connectivity
		% matrices
		CurRegionsI = sub2ind(size(CountA), ToAddStartRegions, ToAddEndRegions);
		
		[CurRegionsISorted, CurRegionsISortedIIDX] = sort(CurRegionsI);
		
		% get the sizes and indices of 
		FirstIndices = [1; find(diff(CurRegionsISorted) > 0) + 1];
		CurSZ = diff([FirstIndices; length(CurRegionsI) + 1]);
		
		LengthAValues = mat2cell_vec(ToAddArcLengths(CurRegionsISortedIIDX), int32(CurSZ));
		WeightedAValues = mat2cell_vec(1 ./ ToAddArcLengths(CurRegionsISortedIIDX), int32(CurSZ));
		
		CurCountA = zeros(NumROIs);
		CurWeightedA = zeros(NumROIs);
		CurLengthA = zeros(NumROIs);
		
		CurCountA(CurRegionsISorted(FirstIndices)) = CurSZ;
		CurWeightedA(CurRegionsISorted(FirstIndices)) = cellfun(@sum, WeightedAValues);
		CurLengthA(CurRegionsISorted(FirstIndices)) = cellfun(@sum, LengthAValues);
		
		CurCountA = CurCountA + CurCountA';
		CurWeightedA = CurWeightedA + CurWeightedA';
		CurLengthA = CurLengthA + CurLengthA';
		
		CountA(:, :, CurIter) = CountA(:, :, CurIter) + CurCountA;
		WeightedA(:, :, CurIter) = WeightedA(:, :, CurIter) + CurWeightedA;
		LengthA(:, :, CurIter) = LengthA(:, :, CurIter) + CurLengthA;
		
		clear CurCountA CurWeightedA CurLengthA WeightedAValues LengthAValues CurSZ FirstIndices CurRegionsISorted CurRegionsISortedIIDX CurRegionsI;
		
		FinalStartRegions(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd, CurIter) = ToAddStartRegions;
		FinalEndRegions(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd, CurIter) = ToAddEndRegions;
        FinalArcLengths(NumTracksSoFar + 1:NumTracksSoFar + NumTracksToAdd, CurIter) = ToAddArcLengths;
        clear ToAddStartRegions ToAddEndRegions ToAddArcLengths;
		NumTracksSoFar = NumTracksSoFar + NumTracksToAdd;
		
		clear CurTracks SeedValues RejectedTracks AcceptedTracks;
		disp(['Number of tracks generated so far: ' num2str(NumTracksSoFar) ' of ' num2str(TotalNumTracks)]);
	end
	
	% not saving any of the tracts, just start tracking again
	
	SeedSizeMatrix = bsxfun(@plus, SeedSizes(:), SeedSizes(:)');
	SeedSizeMatrix = SeedSizeMatrix + double(SeedSizeMatrix == 0);
	SizeWeightedA(:, :, CurIter) = CountA(:, :, CurIter) ./ SeedSizeMatrix;
	WeightedA(:, :, CurIter)  = WeightedA(:, :, CurIter) .* 2 ./ SeedSizeMatrix;
	T = CountA(:, :, CurIter);
	T(T == 0) = 1;
	LengthA(:, :, CurIter) = LengthA(:, :, CurIter) ./ T;
	clear T;
	
end
save(fullfile(Subject, ['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalStartRegions', 'FinalEndRegions', 'WeightedA', 'CountA', 'SizeWeightedA', 'LengthA', 'SeedSizes', 'FinalArcLengths', '-v7.3');
