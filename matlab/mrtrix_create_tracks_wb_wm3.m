function [varargout] = mrtrix_create_tracks_wb_wm3(Subject, MrtrixMethod, SeedType, Curvature)

% test for testing the zero length track, the last in the file

%[CortexLabels] = freesurfer_aparc_labels(['/' fullfile('home', 'addo', 'matlab', 'FreeSurferColorLUT.txt')]);

%[~, CortexLabels] = read_freesurfer_ctx_labels(['/' fullfile('home', 'addo', 'matlab', 'FreeSurferColorLUT.txt')]);
%keyboard;
% do the aparc tracts
%load freesurfer_cortex_labels;
CortexLabels = load_freesurfer_cortex_labels;

% construct a labels image
%WMMaskMIF = load_mif(fullfile(Subject, 'freesurfer_wm.mif'));
%T1NII = load_nii(fullfile('..', 'freesurfer', Subject, 'mri', 'orig', '001'));
%T1IMG = flipdim(permute(T1NII.img, [2 1 3]), 1);
%WMMaskNII = load_nii(fullfile(Subject, 'freesurfer_wm'));
%FANII = load_nii(fullfile(Subject, 'fa'));
%FAIMG = flipdim(permute(FANII.img, [2 1 3]), 1);
%keyboard;
%WMMaskIMG = flipdim(permute(WMMaskMIF.img, [2 1 3]), 1);
%WMMaskIMG = WMMaskMIF.canonicalimg;
%WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);
%WMMaskNII.img = [];

%SavingMIF = WMMaskMIF;
%SavingMIF.layout.directions = {'-', '+', '+'};
%SavingMIF.layout.dimensions = 0:2;
%keyboard;
%WMMaskIMGSZ = int32(size(WMMaskIMG));
%NumInSlice = WMMaskIMGSZ(1) * WMMaskIMGSZ(2);
%ScaledTransform = WMMaskMIF.transform;
%ScaledTransform = WMMaskNII.hdr.transform;
%T = eye(4);
%T(1:3, 1:3) = diag(WMMaskMIF.vox(1:3));
%keyboard;
%ScaledTransform(1:3, 1:3) = ScaledTransform; 0 0 0 1] * T;
%ScaledTransform = ScaledTransform * T;
%InvTransform = inv(ScaledTransform);
% z = 1;
% CurTracks{1} = [-4.56, -22.1, 39.15];
% 
% CurTracks{z} = InvTransform * [CurTracks{1}'; ones(1, size(CurTracks{1}, 1))];
% CurTracks{1} = CurTracks{1}(1:3, :)';
% %CurTracks{1} = CurTracks{1} .* repmat
% CurTracks{1} = CurTracks{1} ./ repmat(WMMaskMIF.vox(1:3)', size(CurTracks{1}, 1), 1);
% % change to 1-indexing
% %CurTracks{1} = CurTracks{z} + 1;
% 
% %CurTracks{1}(:, 1) = double(WMMaskMIF.dim(1)) - CurTracks{z}(:, 1) + 1;
% %CurTracks{1}(:, 2) = double(WMMaskMIF.dim(2)) - CurTracks{z}(:, 2) + 1;
% 
% keyboard;		
%Labels = cell(1, 2);
%Labels{1} = CurLabels.values;
%Labels{2} = CurLabelsOhNine.values;
%Labels{3} = CortexLabels.WMPARC.values;

%SeedType = {'aparc', 'aparc.a2009s', 'wmparc'};
%SeedType = {'aparc.a2009s'};
%SeedType = {'aparc'};

%WMMaskMIF = load_mif(fullfile(Subject, 'freesurfer_wm.mif'));
%WMMaskIMG = WMMaskMIF.canonicalimg;
%WMMaskMIF.img = [];
%WMMaskMIF.canonicalimg = [];

WMMaskNII = load_nii(fullfile(Subject, 'freesurfer_wm'));
%WMMaskIMG = WMMaskNII.img;
WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);
%WMMaskNII.img = [];

% this creates too many tracks for high res volumes
% create two tracks per voxel of 0.75 x 0.75 x 0.75
WMMaskNIIVoxVolume = prod(WMMaskNII.hdr.dime.pixdim(2:4));

DesiredVolume = 0.7 * 0.7 * 0.7;
%V = 0.9 * 0.44 * 0.44;
%keyboard;
TotalNumTracks = round(WMMaskNIIVoxVolume / DesiredVolume * sum(WMMaskIMG(:) > 0));
%disp([Subject ' ' num2str(TotalNumTracks)]);
clear DesiredVolume;
%keyboard;
%TotalNumTracks = sum(WMMaskIMG(:) > 0);

%TotalNumTracks = 50000;
TotalNumTracks = 5000;
NumTracksSoFar = 0;
%FinalTracks = cell(1, TotalNumTracks);
%FinalArcLengths = zeros(1, TotalNumTracks);
FinalStartRegions = zeros(1, TotalNumTracks);
FinalEndRegions = zeros(1, TotalNumTracks);
TracksHeader = [];

%SeedType = 'subdivided1000';
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
	return;
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

%[~, SeedIMG] = ismember(SeedIMG, CurLabels.values);
%[SeedDist, SeedDistNearestIDX] = bwdist(SeedIMG > 0, 'euclidean');
%SeedSizes = histc(SeedIMG(:), 0.5 + (0:length(CurLabels.values)));
SeedSizes = histc(SeedIMG(:), 1:length(CurLabels.values)) * WMMaskNIIVoxVolume;
%keyboard;
if(exist(fullfile(Subject, ['tracks_' MrtrixMethod '.mat']), 'file') == 2)
	TracksAlreadyGenerated = false;
else
	TracksAlreadyGenerated = false;
end

%SameRegionTracks = cell(0);
%RejectedTracks = cell(0);
%AcceptedThenRejectedTracks = cell(0);
FinalTracks = cell(TotalNumTracks, 1);
%FinalOriginalTracks = cell(TotalNumTracks, 1);
%RealStartRegions = zeros(1, TotalNumTracks);
%RealEndRegions = zeros(1, TotalNumTracks);

AcceptedSeedsNII = WMMaskNII;
%AcceptedSeedsMIF.datatype.MATLABType = 'int32';
%RejectedSeedsMIF = SavingMIF;
%RejectedSeedsMIF.datatype.MATLABType = 'int32';
AcceptedSeedsIMG = zeros(size(SeedIMG), 'int32');
%RejectedSeedsIMG = zeros(size(SeedIMG), 'int32');
%FinalAcceptedTracks = cell(0);
%FinalRejectedTracks = cell(0);

%RejectedCodes = [];
%FinalRejectedEndPoints = cell(0);
%RejectedEndsMIF = SavingMIF;
%RejectedEndsIMG = zeros(size(SeedIMG), 'int32');

WeightedA = zeros(length(CurLabels.values));
LengthA = zeros(length(CurLabels.values));
CountA = zeros(length(CurLabels.values));

IncludeRejected = false;

if(IncludeRejected)
	MrtrixDirectory = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9')];
else
	MrtrixDirectory = [filesep fullfile('home', 'addo', 'dev', 'mrtrix-0.2.9-accepted-only')];
end

BlockSize = 50000;
%keyboard;
if(~TracksAlreadyGenerated)
	while(NumTracksSoFar < TotalNumTracks)
	%while(0)
		%Curvature = 1;
		MrtrixStep = 1;
		%TrackArcLengthThresh = 10;
		%NumTracks = 100000; %round((TotalNumTracks - NumTracksSoFar) * 5);
		NumTracks = max(min((TotalNumTracks - NumTracksSoFar) * 2, 2000000), 1000);
		%NumTracks = TotalNumTracks - NumTracksSoFar;
		%NumTracks = 1000;
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
		if(isempty(TracksHeader))
			TracksHeader = MIFTracks;
			TracksHeader.Tracks = [];
		end
		%MIFTracks.Tracks = [];
		
		%MIFTracks.Tracks = tracks_world_to_img(MIFTracks.Tracks, WMMaskNII);%$, SZ);
		%save_trackvis_tracks_nii(fullfile(Subject, 'sdf.trk'), MIFTracks, WMMaskNII);
		%keyboard;
		%%
		%keyboard;
		%for z = 1:size(result, 1)
		%z = 1;
		
%		if(IncludeRejected)
			
			tokens = textscan(result, '%s %f %f %f %d %d %d %d %d');
			%clear result;
			for z = 1:length(tokens)
				tokens{z} = tokens{z}(1:MIFTracks.count);
			end
%		end
		%tokens = cat(1, tokens{:});
		%[tokens] = regexp(result, '(accept|reject): (-?\d+(\.\d+)?) (-?\d+(\.\d+)?) (-?\d+(\.\d+)?) (-?\d+) (\d) (\d) (\d+) (\d+) (\d+) (\d+)', 'tokens');
		%keyboard;
		
		clear MIFTracks;
		%tokens = tokens(1:NumTracks);
		%keyboard;
		%keyboard;
		%tokenNumbers = str2double(tokens(:, 2:end));
		%tokenNumbers = tokens(2:end);
		TracksSZ = cellfun('size', CurTracks, 1);
		CurNumTracks = numel(CurTracks);
		%ArcLengths = cellfun(@arc_length, CurTracks);
		ArcLengths = TracksSZ * MrtrixStep;
		%OriginalTracks = CurTracks;
		
		%CurTracks2 = tracks_world_to_img(CurTracks, WMMaskNII, size(WMMaskIMG));
		% block this to save RAM
		for z = 1:BlockSize:CurNumTracks
			RightIDX = min(z + BlockSize - 1, CurNumTracks);
			CurIDX = (z:RightIDX);
					
			CurTracks(CurIDX) = tracks_world_to_img(CurTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
		end
		% do this in a for loop to save RAM
		% remove any tracks that are outside the volume
		%T = cat(1, CurTracks{:});
		%OutsideIMG = any(T < 1, 2) | T(:, 1) > size(WMMaskIMG, 2) | T(:, 2) > size(WMMaskIMG, 1) | T(:, 3) > size(WMMaskIMG, 3);
		%clear T;
		%sum(OutsideIMG)
		%OutsideIMG = mat2cell_vec(OutsideIMG, TracksSZ);
		%TracksOutsideIMG = cellfun(@any, OutsideIMG);
		
		TracksOutsideIMG = false(size(CurTracks));
		WMMaskIMGSZ = size(WMMaskIMG);
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
% 		TracksOutsideIMG2 = false(size(CurTracks));
% 		for z = 1:length(CurTracks)
% 			TracksOutsideIMG2(z) = any(CurTracks{z}(:) < 1) || ...
% 				any(CurTracks{z}(:, 1) > WMMaskIMGSZ(2)) || ...
% 				any(CurTracks{z}(:, 2) > WMMaskIMGSZ(1)) || ...
% 				any(CurTracks{z}(:, 3) > WMMaskIMGSZ(3));
% 		end
% 		keyboard;
% 		return;
		TracksInsideIMG = ~TracksOutsideIMG;
		if(any(TracksOutsideIMG))
			CurTracks = CurTracks(TracksInsideIMG);
			TracksSZ = TracksSZ(TracksInsideIMG);
			ArcLengths = ArcLengths(TracksInsideIMG);
		end
		
		disp([num2str(sum(TracksOutsideIMG)) ' tracks outside image volume']);
		clear OutsideIMG TracksOutsideIMG;
		%token
		for z = 1:length(tokens)
			tokens{z} = tokens{z}(TracksInsideIMG);
		end
%		if(IncludeRejected)
		AcceptedTrackMask = TracksInsideIMG & strcmp(tokens{1}, 'accept:');
		disp([num2str(sum(AcceptedTrackMask)) ' tracks left after insideimg']);
%		else
%			AcceptedTrackMask = TracksInsideIMG;
%		end
		%CurTracks2 = tracks_world_to_img(MIF2.Tracks, WMMaskNII, size(WMMaskIMG));
		%keyboard;
		clear TracksInsideIMG;
		%AcceptReject = cellfun(@(x) (x{1}), tokens, 'UniformOutput', false);
		%AcceptedTrackMask = strcmp(tokens(:, 1), 'accept');
		
		%clear TracksInsideIMG;
		disp([num2str(sum(AcceptedTrackMask)) ' mrtrix accepted tracks']);
		disp([num2str(sum(~AcceptedTrackMask)) ' mrtrix rejected tracks']);
		%OrigAcceptedTrackMask = strcmp(tokens{1}, 'accept:');
		StartRegion = tokens{8};%;cellfun(@(x) (str2double(x{10})), tokens, 'UniformOutput', true);
		%[~, StartRegion] = ismember(StartRegion, CurLabels.values);
		EndRegion = tokens{9};%cellfun(@(x) (str2double(x{11})), tokens, 'UniformOutput', true);
		%[~, EndRegion] = ismember(EndRegion, CurLabels.values);
		AcceptedTrackMask(StartRegion == EndRegion) = 0;
		disp([num2str(sum(AcceptedTrackMask)) ' tracks left after same region cull']);
		%keyboard;
		NumBeforeNoWM = sum(AcceptedTrackMask);
		
		% block this to save RAM
		%[RejectedTracksNoWM2] = tracks_no_wm_mask(CurTracks(AcceptedTrackMask), WMMaskIMG);%OriginalSeedIMG);
		CurNumTracks = sum(AcceptedTrackMask);
		IDX = find(AcceptedTrackMask);
		RejectedTracksNoWM = false(CurNumTracks, 1);
		for z = 1:BlockSize:CurNumTracks
			RightIDX = min(z + BlockSize - 1, CurNumTracks);
			CurIDX = z:RightIDX;
			RejectedTracksNoWM(CurIDX) = tracks_no_wm_mask(CurTracks(IDX(CurIDX)), WMMaskIMG);%OriginalSeedIMG);
		end
		clear IDX;
		%keyboard;
		AcceptedTrackMask(AcceptedTrackMask) = ~RejectedTracksNoWM;
		clear RejectedTracksNoWM;
		NumBeforeAnteriorPosterior = sum(AcceptedTrackMask);
		disp([num2str(sum(AcceptedTrackMask)) ' tracks left after no_wm']);
		%sum(AcceptedTrackMask)
		% block this to save RAM
		%[RejectedAnteriorCCPosterior] = tracks_anterior_cc_posterior_trajectory(CurTracks(AcceptedTrackMask), OriginalSeedIMG);
		%[RejectedAnteriorCCPosterior2] = tracks_anterior_cc_posterior_trajectory(CurTracks(AcceptedTrackMask), OriginalSeedIMG);
		CurNumTracks = sum(AcceptedTrackMask);
		IDX = find(AcceptedTrackMask);
		
		RejectedAnteriorCCPosterior = false(CurNumTracks, 1);
		for z = 1:BlockSize:CurNumTracks
			RightIDX = min(z + BlockSize - 1, CurNumTracks);
			CurIDX = z:RightIDX;
			RejectedAnteriorCCPosterior(CurIDX) = tracks_anterior_cc_posterior_trajectory(CurTracks(IDX(CurIDX)), OriginalSeedIMG);
		end
		AcceptedTrackMask(AcceptedTrackMask) = ~RejectedAnteriorCCPosterior;
		%keyboard;
		clear RejectedAnteriorCCPosterior;
		%[RejectedAnteriorCCPosterior] = tracks_anterior_cc_posterior_trajectory(RejectedTracks, OriginalSeedIMG);
		disp([num2str(NumBeforeNoWM - NumBeforeAnteriorPosterior) ' culled by no wm']);
		disp([num2str(NumBeforeAnteriorPosterior - sum(AcceptedTrackMask)) ' culled by anterior posterior']);
		AcceptedTracks = CurTracks(AcceptedTrackMask);
%		AcceptedOriginalTracks = OriginalTracks(AcceptedTrackMask);
		%RejectedTracks = CurTracks(~AcceptedTrackMask);
		%FinalRejectedTracks = [FinalRejectedTracks; CurTracks(~AcceptedTrackMask)];
		%keyboard;
		%AcceptedArcLengths = ArcLengths(AcceptedTrackMask);
		%AcceptedStartRegions = StartRegion(AcceptedTrackMask);
		%AcceptedEndRegions = EndRegion(AcceptedTrackMask);
		
		%clear F;
		%$RejectedTracksIDX = find(~AcceptedTrackMask);
		
		%AcceptedSeeds = cellfun(@(x) (cat(2, str2double(x{2}),  str2double(x{3}),  str2double(x{4}))), tokens(AcceptedTrackMask), 'UniformOutput', false);
		%AcceptedSeeds = cat(1, AcceptedSeeds{:});
		Seeds = cat(2, tokens{2:4});
		AcceptedSeeds = Seeds(AcceptedTrackMask, :);
		%RejectedSeeds = Seeds(~AcceptedTrackMask, :);
		clear Seeds;
% 		RejectedSeeds = cellfun(@(x) (cat(2, str2double(x{2}),  str2double(x{3}),  str2double(x{4}))), tokens(~AcceptedTrackMask), 'UniformOutput', false);
% 		RejectedSeeds = cat(1, RejectedSeeds{:});
		%RejectedSeeds = tokenNumbers(~AcceptedTrackMask, 1:3);
		%Codes = cellfun(@(x) (str2double(x{5})), tokens(~AcceptedTrackMask), 'UniformOutput', true);
		Codes = tokens{5};
		%RejectedCodes = [RejectedCodes; Codes(~AcceptedTrackMask)];
		%CurRejectedCodes = Codes(~AcceptedTrackMask);
		UniqueCodes = unique(Codes);
		for z = 1:length(UniqueCodes)
			disp([num2str(sum(Codes == UniqueCodes(z))) ' got code ' num2str(UniqueCodes(z))]);
		end

		if(IncludeRejected)
			Excluded = logical(tokens{6});
			Included = logical(tokens{7});
			disp([num2str(sum(~Excluded)) ' exclusion region']);
			disp([num2str(sum(~Included)) ' didnt reach include regions']);
			disp([num2str(sum(TckSize < MinSize)) ' didnt reach required size']);
		end
		clear tokens;
		%Excluded = logical(tokenNumbers(~AcceptedTrackMask, 5));
		%Included = logical(tokenNumbers(~AcceptedTrackMask, 6));
		%TckSize = tokenNumbers(~AcceptedTrackMask, 7);
		%MinSize = tokenNumbers(~AcceptedTrackMask, 8);
		
% 		Excluded = cellfun(@(x) (str2double(x{6})), tokens(~AcceptedTrackMask), 'UniformOutput', true);
% 		Included = cellfun(@(x) (str2double(x{7})), tokens(~AcceptedTrackMask), 'UniformOutput', true);
% 		TckSize = cellfun(@(x) (str2double(x{8})), tokens(~AcceptedTrackMask), 'UniformOutput', true);
% 		MinSize = cellfun(@(x) (str2double(x{9})), tokens(~AcceptedTrackMask), 'UniformOutput', true);
		
% 		StartRegion = cellfun(@(x) (str2double(x{10})), tokens, 'UniformOutput', true);
% 		EndRegion = cellfun(@(x) (str2double(x{11})), tokens, 'UniformOutput', true);
		
		disp([num2str(sum(~AcceptedTrackMask)) ' rejected tracks']);
		
		%keyboard;
%		disp([num2str(sum(Codes == -4)) ' rejected due to curvature']);
% 		disp([num2str(sum(Codes == -4 & TckSize < MinSize)) ' didnt reach size due to curvature']);
% 		disp([num2str(sum(Codes == -4 & TckSize > MinSize & ~Included)) ' reached size but didnt hit includes due to curvature']);
% 		
		%%
		
		%CurTracks = load_strands(fullfile(Subject, ['tracks_' CurSeedType], ['tracks_' num2str(CurLabels(CurLabel)) '.tck.gz']));
		%CurTracks = load_strands(TrackFile);

		%length(CurTracks)
		%length(tokens)
		%keyboard;
		%AcceptedTracks = CurTracks(AcceptedTrackMask);
		%RejectedTracks = CurTracks(~AcceptedTrackMask);
		%CurTracks = CurTracks(AcceptedTracks);
		%[FornixTracks, AnyInterAccepted, AnyInterNoCCAccepted, AnyInterAnteriorCCAccepted, AnyInterPosteriorCCAccepted] = tracks_through_fornix(AcceptedTracks, OriginalSeedIMG);
		%[FornixTracks, AnyInterRejected, AnyInterNoCCRejected, AnyInterAnteriorCCRejected, AnyInterPosteriorCCRejected] = tracks_through_fornix(RejectedTracks, OriginalSeedIMG);
		
		%keyboard;
		%clear MIFTracks;
		
		%if(strcmp(Subject, 'IRC143-1001'))
		%	copyfile(TrackFile, fullfile(Subject, 'tmp4.tck'));
		%end
		%keyboard;
		%delete(TrackFile);
		
		%TracksSZ = cellfun(@(x) (size(x, 1)), AcceptedTracks);
		%%CurTracks = cat(1, CurTracks{:});
		%length(CurTracks)
		
		%CumulativeArcLengths = cell(size(CurTracks));
		%disp(num2str(length(CurTracks)));
		%TracksSZ
		%CurTracks = tracks_world_to_img(CurTracks, WMMaskMIF, size(WMMaskIMG));
		%AcceptedTracks = tracks_world_to_img(AcceptedTracks, WMMaskMIF, size(WMMaskIMG));
		%RejectedTracks = tracks_world_to_img(RejectedTracks, WMMaskMIF, size(WMMaskIMG));
		%CurRejectedCodes = Codes(~AcceptedTrackMask);
		%CurRejectedTckSize = TckSize(~AcceptedTrackMask);
		%CurRejectedMinSize = MinSize(~AcceptedTrackMask);
% 		DoNearlyAccepted = 0;
% 		if(DoNearlyAccepted && IncludeRejected)
% 
% 			RejectedStartEndCoords = cell(length(RejectedTracks), 1);
% 			StartVectors = cell(length(RejectedTracks), 1);
% 			EndVectors = cell(length(RejectedTracks), 1);
% 
% 			RejectedTracksSZ = TracksSZ(~AcceptedTrackMask);
% 
% 			for z = 1:length(RejectedTracks)
% 				RejectedStartEndCoords{z} = RejectedTracks{z}([1 RejectedTracksSZ(z)], :);
% 				if(RejectedTracksSZ(z) > 2)
% 					StartVectors{z} = RejectedTracks{z}(1, :) - RejectedTracks{z}(2, :);
% 					EndVectors{z} = RejectedTracks{z}(RejectedTracksSZ(z), :) - RejectedTracks{z}(RejectedTracksSZ(z) - 1, :);
% 				else
% 					StartVectors{z} = [0 0 0];
% 					EndVectors{z} = [0 0 0];
% 				end
% 			end
% 
% 			%[RejectedStartEndCoords] = cellfun(@first_last_coords, RejectedTracks, 'UniformOutput', false);
% 			RejectedStartEndCoords = cat(3, RejectedStartEndCoords{:});
% 			%keyboard;
% 			%M = CurRejectedCodes < 0 & CurRejectedTckSize > CurRejectedMinSize;
% 		
% 			CurTrackSeedDists = squeeze(interp3_linear_fast(SeedDist, RejectedStartEndCoords(:, 1, :), RejectedStartEndCoords(:, 2, :), RejectedStartEndCoords(:, 3, :), 'linear'));
% 			NearestPoints = squeeze(interp3_linear_fast(SeedDistNearestIDX, RejectedStartEndCoords(:, 1, :), RejectedStartEndCoords(:, 2, :), RejectedStartEndCoords(:, 3, :), 'nearest'));
% 			%[StartVectors, EndVectors] = cellfun(@streamline_start_end_vectors, RejectedTracks, 'UniformOutput', false);
% 			StartVectors = cat(1, StartVectors{:});
% 			StartVectorsMAG = sqrt(sum(StartVectors .* StartVectors, 2));
% 			StartVectorsMAG(StartVectorsMAG == 0) = 1;
% 			StartVectors = bsxfun(@rdivide, StartVectors, StartVectorsMAG);
% 
% 			EndVectors = cat(1, EndVectors{:});
% 			EndVectorsMAG = sqrt(sum(EndVectors .* EndVectors, 2));
% 			EndVectorsMAG(EndVectorsMAG == 0) = 1;
% 			EndVectors = bsxfun(@rdivide, EndVectors, EndVectorsMAG);
% 			clear StartVectorsMAG EndVectorsMAG;
% 
% 			[NearestPointsI, NearestPointsJ, NearestPointsK] = ind2sub(size(SeedIMG), NearestPoints);
% 			NearestStartPoints = [NearestPointsJ(1, :)', NearestPointsI(1, :)', NearestPointsK(1, :)'];
% 			NearestEndPoints = [NearestPointsJ(2, :)', NearestPointsI(2, :)', NearestPointsK(2, :)'];
% 			clear NearestPointsI NearestPointsJ NearestPointsK;
% 			%keyboard;
% 			%[StartOfStreamlines, EndOfStreamlines] = cellfun(@streamline_start_end_points, RejectedTracks, 'UniformOutput', false);
% 			%StartOfStreamlines = cellfun(@(x) (x(1, :)), RejectedTracks, 'UniformOutput', false);
% 			%StartOfStreamlines = cat(1, StartOfStreamlines{:});
% 			%EndOfStreamlines = cellfun(@(x) (x(end, :)), RejectedTracks, 'UniformOutput', false);
% 			%EndOfStreamlines = cat(1, EndOfStreamlines{:});
% 
% 			StartToNearest = NearestStartPoints - squeeze(RejectedStartEndCoords(1, :, :))';
% 			EndToNearest = NearestEndPoints - squeeze(RejectedStartEndCoords(2, :, :))';
% 
% 			StartToNearestMAG = sqrt(sum(StartToNearest .* StartToNearest, 2));
% 			EndToNearestMAG = sqrt(sum(EndToNearest .* EndToNearest, 2));
% 			StartToNearestMAG(StartToNearestMAG == 0) = 1;
% 			EndToNearestMAG(EndToNearestMAG == 0) = 1;
% 			StartToNearest = bsxfun(@rdivide, StartToNearest, StartToNearestMAG);
% 			EndToNearest = bsxfun(@rdivide, EndToNearest, EndToNearestMAG);
% 			clear StartToNearestMAG EndToNearestMAG;
% 
% 			StartDotProducts = sum(StartToNearest .* StartVectors, 2);
% 			EndDotProducts = sum(EndToNearest .* EndVectors, 2);
% 
% 			DegTHRESH = cos(90 / 180 * pi);
% 			%keyboard;
% 			%any(isnan(NearestPoints(:)))
% 			StartRegion(~AcceptedTrackMask) = SeedIMG(NearestPoints(1, :));
% 			EndRegion(~AcceptedTrackMask) = SeedIMG(NearestPoints(2, :));
% 
% 			M = CurRejectedCodes ~= 8 & ...
% 				StartDotProducts > DegTHRESH & EndDotProducts > DegTHRESH & ...
% 				all(CurTrackSeedDists < 2)' & ...
% 				StartRegion(~AcceptedTrackMask) ~= EndRegion(~AcceptedTrackMask) & ...
% 				StartRegion(~AcceptedTrackMask) > 0 & ...
% 				EndRegion(~AcceptedTrackMask) > 0;
% 			clear DotProductsStart DotProductsEnd EndToNearest StartToNearest EndOfStreamlines StartOfStreamlines;
% 			clear NearestEndPoints NearestStartPoints StartVectors EndVectors CurTrackSeedDists NearestPoints CurRejectedCodes;
% 
% 			AcceptedTracks = [AcceptedTracks; RejectedTracks(M)];
% 			RejectedTracks = RejectedTracks(~M);
% 			AcceptedSeeds = [AcceptedSeeds; RejectedSeeds(M, :)];
% 			RejectedSeeds = RejectedSeeds(~M);
% 
% 			AcceptedTrackMask(~AcceptedTrackMask) = M;
% 			clear M;
% 		end
		
		AcceptedArcLengths = ArcLengths(AcceptedTrackMask);
		AcceptedStartRegions = StartRegion(AcceptedTrackMask);
		AcceptedEndRegions = EndRegion(AcceptedTrackMask);
		clear RejectedTracks;
		
 		%keyboard;
% 		for z = 1:length(RejectedTracks)
% 			if(CurRejectedCodes(z) < 0 && CurRejectedTckSize(z) > CurRejectedMinSize(z))
% 				%%
% 				CurTrackSeedDists = interp3(SeedDist, RejectedTracks{z}([1 end], 1), RejectedTracks{z}([1 end], 2), RejectedTracks{z}([1 end], 3), 'linear');
% 				NearestPoints = interp3(SeedDistNearestIDX, RejectedTracks{z}([1 end], 1), RejectedTracks{z}([1 end], 2), RejectedTracks{z}([1 end], 3), 'nearest');
% 				[ID, JD, KD] = ind2sub(size(SeedIMG), NearestPoints);
% 				StartVector = RejectedTracks{z}(1, :) - RejectedTracks{z}(2, :);
% 				StartVector = StartVector ./ sqrt(sum(StartVector .* StartVector));
% 				EndVector = RejectedTracks{z}(end, :) - RejectedTracks{z}(end - 1, :);
% 				EndVector = EndVector ./ sqrt(sum(EndVector .* EndVector));
% 				StartToNearestVector = [JD(1) - RejectedTracks{z}(1, 1), ID(1) - RejectedTracks{z}(1, 2), KD(1) - RejectedTracks{z}(1, 3)];
% 				StartToNearestVector = StartToNearestVector ./ sqrt(sum(StartToNearestVector .* StartToNearestVector));
% 				EndToNearestVector = [JD(2) - RejectedTracks{z}(end, 1), ID(2) - RejectedTracks{z}(end, 2), KD(2) - RejectedTracks{z}(end, 3)];
% 				EndToNearestVector = EndToNearestVector ./ sqrt(sum(EndToNearestVector .* EndToNearestVector));
% 				
% 				
% 				DotProductStart = sum(StartVector .* StartToNearestVector);
% 				DotProductEnd = sum(EndVector .* EndToNearestVector);
% 				%%
% 				keyboard;
% 			end
% 		end
% 		
%		for z = 1:length(RejectedTracks)
			%D = diff(AcceptedTracks{z}, 1, 1);
			%CumulativeArcLength{z} = [0; cumsum(sqrt(sum(D .* D, 2)))];
			%clear D;
% 			RejectedTracks{z} = InvTransform * [RejectedTracks{z}'; ones(1, size(RejectedTracks{z}, 1))];
% 			RejectedTracks{z} = RejectedTracks{z}(1:3, :)';
% 			RejectedTracks{z} = bsxfun(@rdivide, RejectedTracks{z}, WMMaskMIF.vox(1:3)');
% 			% change to 1-indexing
% 			RejectedTracks{z} = RejectedTracks{z} + 1;
% 
% 			RejectedTracks{z}(:, 1) = double(WMMaskMIF.dim(1)) - RejectedTracks{z}(:, 1) + 1;
% 			RejectedTracks{z}(:, 2) = double(WMMaskMIF.dim(2)) - RejectedTracks{z}(:, 2) + 1;
% 			
%			I = round(RejectedTracks{z}([1 end], 2)) + size(RejectedEndsIMG, 1) * round(RejectedTracks{z}([1 end], 1)) + size(RejectedEndsIMG, 1) * size(RejectedEndsIMG, 2) * round(RejectedTracks{z}([1 end], 3));
%			for k = 1:length(I)
%				RejectedEndsIMG(I(k)) = RejectedEndsIMG(I(k)) + 1;
%			end
%		end
		%keyboard;
		%AcceptedSeeds = tracks_world_to_img({AcceptedSeeds}, WMMaskMIF, size(WMMaskIMG));
		%AcceptedSeeds = AcceptedSeeds{1};
		
		I = round(AcceptedSeeds(:, 2)) + size(AcceptedSeedsIMG, 1) * round(AcceptedSeeds(:, 1)) + size(AcceptedSeedsIMG, 1) * size(AcceptedSeedsIMG, 2) * round(AcceptedSeeds(:, 3));
		AddToIMG = histc(I, 1:numel(AcceptedSeedsIMG));
		AcceptedSeedsIMG = AcceptedSeedsIMG + int32(reshape(AddToIMG, size(AcceptedSeedsIMG)));
		clear I AddToIMG;
		%for z = 1:length(I)
		%	AcceptedSeedsIMG(I) = AcceptedSeedsIMG(I) + 1;
		%end
%		RejectedSeeds = tracks_world_to_img({RejectedSeeds}, WMMaskMIF, size(WMMaskIMG));
%		RejectedSeeds = RejectedSeeds{1};
		
%		I = round(RejectedSeeds(:, 2)) + size(RejectedSeedsIMG, 1) * round(RejectedSeeds(:, 1)) + size(RejectedSeedsIMG, 1) * size(RejectedSeedsIMG, 2) * round(RejectedSeeds(:, 3));
%		AddToIMG = histc(I, 1:numel(RejectedSeedsIMG));
%		RejectedSeedsIMG = RejectedSeedsIMG + int32(reshape(AddToIMG, size(RejectedSeedsIMG)));
%		clear I AddToIMG;
		
		% all tracks are now in IMAGE coordinates
		
		% go through the accepted tracts and mark the ones that go straight
		% from left WM to right WM or vice versa as rejected, these are
		% tracts going through the fornix
% 		AcceptedStartRegion = StartRegion(AcceptedTrackMask);
% 		AcceptedEndRegion = EndRegion(AcceptedTrackMask);
% % 		AcceptedInterHemispheric = (AcceptedStartRegion < 2000 & AcceptedEndRegion > 2000) | ...
% % 			(AcceptedStartRegion > 2000 & AcceptedEndRegion < 2000);
% 		
% 		RejectedStartRegion = StartRegion(~AcceptedTrackMask);
% 		RejectedEndRegion = EndRegion(~AcceptedTrackMask);
% 		
% 		RejectedInterHemispheric = (RejectedStartRegion < 2000 & RejectedEndRegion > 2000) | ...
% 			(RejectedStartRegion > 2000 & RejectedEndRegion < 2000);
		
		%[FornixTracks, InterhemisphericNoCC, AnyInterAccepted, AnyInterNoCCAccepted] = tracks_through_fornix(AcceptedTracks, OriginalSeedIMG, AcceptedInterHemispheric);
		%[RejectedFornixTracks, RejectedInterhemisphericNoCC, AnyInter, AnyInterNoCC] = tracks_through_fornix(RejectedTracks, OriginalSeedIMG, true(size(RejectedInterHemispheric)));
		%[FornixTracks, InterhemisphericNoCC] = tracks_through_fornix(AcceptedTracks, OriginalSeedIMG, AcceptedInterHemispheric);
		%InterhemisphericNoCC = tracks_interhemispheric_no_cc(AcceptedTracks(AcceptedInterHemispheric), OriginalSeedIMG);

% 		FornixTracksIntraHemispheric = FornixTracks(:) & (...
% 			(AcceptedStartRegion(:) < 2000 & AcceptedEndRegion(:) < 2000) | ...
% 			(AcceptedStartRegion(:) > 2000 & AcceptedEndRegion(:) > 2000));
% 		
% 		
		%A = strfind(CortexLabels.APARC.labels, 'temporal');
		
		%TemporalLabels = ~cellfun(@isempty, A);
		%TemporalValues = CortexLabels.APARC.values(TemporalLabels);
				
		
		%AcceptedTracks = AcceptedTracks(~AcceptedAnteriorCCPosterior);
		% more efficient way of doing this, avoid looping
		% need to make cell arrays of the values at each 
% 		for CurTrackIDX = 1:length(AcceptedTracks)
% 			CountA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) = CountA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) + 1;
% 			CountA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) = CountA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) + 1;
% 			LengthA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) = LengthA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) + AcceptedArcLengths(CurTrackIDX);
% 			LengthA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) = LengthA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) + AcceptedArcLengths(CurTrackIDX);
% 			WeightedA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) = WeightedA(AcceptedStartRegions(CurTrackIDX), AcceptedEndRegions(CurTrackIDX)) + 1 ./ AcceptedArcLengths(CurTrackIDX);
% 			WeightedA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) = WeightedA(AcceptedEndRegions(CurTrackIDX), AcceptedStartRegions(CurTrackIDX)) + 1 ./ AcceptedArcLengths(CurTrackIDX);
% 		end
% 		
		% get the linear indices of the regions in the connectivity
		% matrices
		CurRegionsI = sub2ind(size(CountA), AcceptedStartRegions, AcceptedEndRegions);
		%CountAValues = ones(length(AcceptedTracks), 1);
		%LengthAValues = AcceptedArcLengths;
		%WeightedAValues = 1 ./ AcceptedArcLengths;
		% sort the indices so we know which elements to combine
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
		%%
% 		clf;
% 		%slice(double(ismember(OriginalSeedIMG, TemporalValues)), [], [], 114); axis equal ij;
% 		V = double(ismember(OriginalSeedIMG, CortexLabels.APARC.values));
% 		V = double(ismember(OriginalSeedIMG, TemporalValues));
% 		V(ismember(OriginalSeedIMG, [2 41])) = 2;
% 		h = slice(V, [], [], 120); axis equal ij off; colormap gray;
% 		set(h, 'EdgeColor', 'none');
% 		%streamline(AcceptedTracks(InterhemisphericNoCC))
% 		h = streamline(RejectedTracks(RejectedFornixTracks));
% 		set(h, 'Color', 'r');
% 		h = streamline(RejectedTracks(AnyInter & ~Excluded));
% 		set(h, 'Color', 'b');
% 		h = streamline(AcceptedTracks(AnyInterAccepted));
% 		set(h, 'Color', 'g');
% 		view(-6, 62);
% 		view(-52, 20);
		%h = streamline(AcceptedTracks);
		%set(h, 'Color', 'g');
		%%
		%for z = 1:length(I)
		%	RejectedSeedsIMG(I) = RejectedSeedsIMG(I) + 1;
		%end
%		keyboard;
		%AcceptedTracks = mat2cell(AcceptedTracks, TracksSZ, 3);
		
		%%
		%StartRegions = zeros(1, length(CurTracks));
		%EndRegions = zeros(1, length(CurTracks));
		%ValidTracks = false(1, length(CurTracks));

		%TracksSZ = cellfun(@numrows, CurTracks);

		%T = cat(1, AcceptedTracks{:});
		
		%SeedValues = interp3_nearest_non_zero(int32(SeedIMG), double(T(:, 2)), double(T(:, 1)), double(T(:, 3)));
		%SeedValues = interp3(double(OriginalSeedIMG), double(T(:, 2)), double(T(:, 1)), double(T(:, 3)), 'nearest');
		%SeedValues = mat2cell_vec(SeedValues, TracksSZ);
		%keyboard;
		%XYZInIMG = all(T >= 1, 2) & T(:, 1) <= WMMaskIMGSZ(2) & T(:, 2) <= WMMaskIMGSZ(1) & T(:, 3) <= WMMaskIMGSZ(3);
		%CurValidTracks = 1:length(CurTracks);
		%XYZInIMG = mat2cell(XYZInIMG, TracksSZ, 1);
		%XYZInIMG = mat2cell_vec(XYZInIMG, TracksSZ);
		%TracksInIMG = cellfun(@all, XYZInIMG);
		%keyboard;
		%TracksInIMG = TracksInIMG & (ArcLengths > 10);
		%clear XYZInIMG;

		%CurValidTracks = TracksInIMG;
		%ValidTrackIDX = find(CurValidTracks);
		
		%clear T;
		
		%CurRejectedTracks = false(size(CurTracks));
		%CurSameRegionTracks = false(size(CurTracks));
		%StartAndEnd = zeros(2, length(TracksSZ));
		%CurAcceptedTracks = 0;
		%for z = 1:length(CurTracks)
			% new strategy, look for non-zero values separated by runs of
			% zeros, compute streamlines from those
			
% 			F = find(SeedValues{z} > 0);
% 			if(~isempty(F))
% 				for k = 1:length(F) - 1
% 					%keyboard;
% 					if((F(k + 1) - F(k)) * MrtrixStep > TrackArcLengthThresh && SeedValues{z}(F(k)) ~= SeedValues{z}(F(k + 1)))
% 						FinalTracks{NumTracksSoFar + 1} = CurTracks{z}(F(k):F(k + 1), :);
% 						FinalStartRegions(NumTracksSoFar + 1) = SeedValues{z}(F(k));
% 						FinalEndRegions(NumTracksSoFar + 1) = SeedValues{z}(F(k + 1));
% 						NumTracksSoFar = NumTracksSoFar + 1;
% 					end
% 				end
%  			end
			%keyboard;
			%NonZeroEntries = SeedValues{z} > 0; 
			%FinalTracks
% 			if(any(SeedValues{z}))
% 				[FirstI] = find(SeedValues{z} > 0, 1, 'first');
% 				[LastI] = find(SeedValues{z} > 0, 1, 'last');
% 
% 				if(FirstI < 5 && LastI > length(SeedValues{z}) - 5)% && SeedValues{z}(FirstI) ~= SeedValues{z}(LastI))
% 					FinalTracks{NumTracksSoFar + 1} = CurTracks{z};
% 					FinalStartRegions(NumTracksSoFar + 1) = SeedValues{z}(FirstI);
% 					FinalEndRegions(NumTracksSoFar + 1) = SeedValues{z}(LastI);
% 					NumTracksSoFar = NumTracksSoFar + 1;
% 					if(NumTracksSoFar == TotalNumTracks)
% 						break;
% 					end
% 				else
% 					% use nearest labels to the start and end points
% 					StartI = round(CurTracks{z}(1, 2)) + ...
% 						size(WMMaskIMG, 1) * round(CurTracks{z}(1, 1)) + ...
% 						size(WMMaskIMG, 1) * size(WMMaskIMG, 2) * round(CurTracks{z}(1, 3));
% 					EndI = round(CurTracks{z}(end, 2)) + ...
% 						size(WMMaskIMG, 1) * round(CurTracks{z}(end, 1)) + ...
% 						size(WMMaskIMG, 1) * size(WMMaskIMG, 2) * round(CurTracks{z}(end, 3));
% 					if(SeedIMG(StartI) == 0)
% 						
% 					end
% 					if(SeedIMG(EndI) == 0)
% 						
% 					end
% 					
% 					%keyboard;
% 					%AcceptedThenRejectedTracks = [AcceptedThenRejectedTracks; CurTracks{z}];
% 				end
% 			
% 			end
% 			[FirstI] = find(SeedValues{z} > 0, 1, 'first');
% 			[LastI] = find(SeedValues{z} > 0, 1, 'last');
% 			if(~isempty(FirstI) && ~isempty(LastI))
% 				if(FirstI ~= LastI && SeedValues{z}(FirstI) ~= SeedValues{z}(LastI))
% 					StartAndEnd(:, z) = [SeedValues{z}(FirstI); SeedValues{z}(LastI)];
% 					CurAcceptedTracks = CurAcceptedTracks + 1;
% 				else
% 					CurSameRegionTracks(z) = 1;
% 				end
% 			else
% 				CurRejectedTracks(z) = 1;
% 			end
		%end
		%keyboard;
		%clear T;
		%TracksInIMGIDX = find(TracksInIMG);
		% WMInterp = interp3(double(WMMaskIMG), T(:, 2), T(:, 1), T(:, 3), 'linear');
		% clear T;
		% WMInterp = (WMInterp >= 0.25);
		% WMInterp = mat2cell(WMInterp(:), TracksSZ, 1);
		
		%[ID, JD, KD] = cellfun(@first_last_coords_to_idx, CurTracks(TracksInIMG), 'UniformOutput', false);
		
		%ID = cat(2, ID{:});
		%JD = cat(2, JD{:});
		%KD = cat(2, KD{:});
		
		%ID = double(ID);
		%JD = double(JD);
		%KD = double(KD);
		
		%StartAndEnd = interp3_nearest_non_zero(int32(SeedIMG), ID, JD, KD);
		%StartAndEnd = reshape(StartAndEnd, size(ID));
		%keyboard;
		%clear ID JD KD;
		%ValidStartAndEndIDX = TracksInIMGIDX(all(StartAndEnd > 0));
		
		%CurValidTracks(any(StartAndEnd == 0)) = 0;
		%StartRegions = StartAndEnd(1, all(StartAndEnd > 0));
		%EndRegions = StartAndEnd(2, all(StartAndEnd > 0));
		
		%ValidTrackIDX = find(CurValidTracks);
		%clear StartAndEnd;
		%for CurTrackIDX = 1:size(StartRegions, 2)
		%	if(StartRegions(CurTrackIDX) ~= EndRegions(CurTrackIDX))
% 				if(StartRegions(CurTrackIDX) > 66 || EndRegions(CurTrackIDX) > 66)
% 					keyboard;
% 				end
		%		RealStartRegions(ValidTrackIDX(CurTrackIDX)) = StartRegions(CurTrackIDX);
		%		RealEndRegions(ValidTrackIDX(CurTrackIDX)) = EndRegions(CurTrackIDX);
		%	else
		%		CurValidTracks(ValidTrackIDX(CurTrackIDX)) = 0;
		%	end
		%end
% 		for CurTrackIDX = find(TracksInIMG)'%1:length(CurTracks)
% 			CurTrackSZ = TracksSZ(CurTrackIDX);%size(CurTracks{CurTrackIDX}, 1);
% 		 	if(all(CurTracks{CurTrackIDX}(:) >= 1) && ...
% 		 			all(CurTracks{CurTrackIDX}(:, 1) < size(SeedIMG, 2)) && ...
% 		 			all(CurTracks{CurTrackIDX}(:, 2) < size(SeedIMG, 1)) && ...
% 		 			all(CurTracks{CurTrackIDX}(:, 3) < size(SeedIMG, 3)))
% 
% 			
% 			ID = int32(floor(CurTracks{CurTrackIDX}([1, CurTrackSZ], 2)));
% 			JD = int32(floor(CurTracks{CurTrackIDX}([1, CurTrackSZ], 1)));
% 			KD = int32(floor(CurTracks{CurTrackIDX}([1, CurTrackSZ], 3)));
% 
% 			ID(ID == size(SeedIMG, 1)) = ID(ID == size(SeedIMG, 1)) - 1;
% 			JD(JD == size(SeedIMG, 2)) = JD(JD == size(SeedIMG, 2)) - 1;
% 			KD(KD == size(SeedIMG, 3)) = KD(KD == size(SeedIMG, 3)) - 1;
% 
% 			%I = sub2ind(size(SeedIMG), ID, JD, KD);
% 			% inline sub2ind
% 
% 			IDValues = [ID(:)'    ; ID(:)' + 1; ID(:)'    ; ID(:)' + 1; ID(:)'    ; ID(:)' + 1; ID(:)'    ; ID(:)' + 1];
% 			JDValues = [JD(:)'    ; JD(:)'    ; JD(:)' + 1; JD(:)' + 1; JD(:)'    ; JD(:)'    ; JD(:)' + 1; JD(:)' + 1];
% 			KDValues = [KD(:)'    ; KD(:)'    ; KD(:)'    ; KD(:)'    ; KD(:)' + 1; KD(:)' + 1; KD(:)' + 1; KD(:)' + 1];
% 
% 			I = IDValues + WMMaskIMGSZ(1) * (JDValues - 1) + NumInSlice * (KDValues - 1);
% 			%StartRegions(CurTrackIDX) = SeedIMG(I(:, 1));
% 			%EndRegions(CurTrackIDX) = SeedIMG(I(:, 2));
% 
% 			CurStartRegions = SeedIMG(I(:, 1));
% 			CurEndRegions = SeedIMG(I(:, 2));
% 
% 			StartTF = (CurStartRegions > 0);
% 			EndTF = (CurEndRegions > 0);
% 			if(any(StartTF) && any(EndTF))
% 
% 				F = CurStartRegions(StartTF);
% 				StartRegions(CurTrackIDX) = F(1);
% 				F = CurEndRegions(EndTF);
% 				EndRegions(CurTrackIDX) = F(1);
% 				clear F;
% 				if(StartRegions(CurTrackIDX) ~= EndRegions(CurTrackIDX))
% 					CountA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) = CountA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) + 1;
% 					CountA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) = CountA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) + 1;
% 					LengthA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) = LengthA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) + ArcLengths(CurTrackIDX);
% 					LengthA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) = LengthA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) + ArcLengths(CurTrackIDX);
% 					WeightedA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) = WeightedA(StartRegions(CurTrackIDX), EndRegions(CurTrackIDX)) + 1 ./ ArcLengths(CurTrackIDX);
% 					WeightedA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) = WeightedA(EndRegions(CurTrackIDX), StartRegions(CurTrackIDX)) + 1 ./ ArcLengths(CurTrackIDX);
% 				else
% 					CurValidTracks(CurTrackIDX) = 0;
% 				end
% 			else
% 				CurValidTracks(CurTrackIDX) = 0;
% 			end
% 
% 			%CurRegions = SeedIMG(I);
% 			%CurRegions = CurRegions(CurRegions > 0);
% 
% 		% 	end
% 
% 		end
		%clear I ID JD KD IDValues JDValues KDValues;

		%FinalRejectedTracks = [FinalRejectedTracks; RejectedTracks];
		%FinalAcceptedTracks = [FinalAcceptedTracks; CurTracks(AcceptedTrackMask)];
		
		TracksToAdd = 1:min(length(AcceptedTracks), TotalNumTracks - NumTracksSoFar);%, length(CurTracks));
		
		%FinalOriginalTracks(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedOriginalTracks(TracksToAdd);
		FinalTracks(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedTracks(TracksToAdd);
		FinalArcLengths(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedArcLengths(TracksToAdd);
		
		FinalStartRegions(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedStartRegions(TracksToAdd);
		FinalEndRegions(NumTracksSoFar + 1:NumTracksSoFar + length(TracksToAdd)) = AcceptedEndRegions(TracksToAdd);
		NumTracksSoFar = NumTracksSoFar + length(TracksToAdd);
		clear CurTracks SeedValues RejectedTracks AcceptedTracks;
		disp(['Number of tracks generated so far: ' num2str(NumTracksSoFar) ' of ' num2str(TotalNumTracks)]);
		%RejectedTracks = [RejectedTracks; CurTracks(CurRejectedTracks)];
		%SameRegionTracks = [SameRegionTracks; CurTracks(CurSameRegionTracks)];
		%ValidTracks(StartRegions == EndRegions) = 0;
		%MIFTracks = load_mif_tracks(fullfile(Subject, 'tracks_wb.tck'));
	end
	
	%load(fullfile(Subject, ['tracks_wb_eddy_cortex_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature)]), '-mat', 'FinalTracks', 'FinalArcLengths', 'FinalStartRegions', 'FinalEndRegions');
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
%	keyboard;
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
	
	%TracksHeader.Tracks = tracks_world_to_img(TracksHeader.Tracks, FANII, SZ);
	%save_trackvis_tracks_nii([F '.trk'], TracksHeader, FANII);
	TracksHeader = EmptyTracksHeader;
	
	% FinalTracks go to trackvis img space
	%FinalTracks = tracks_world_to_img_trackvis(FinalTracks, WMMaskNII, size(WMMaskIMG));
	
% 	for z = 1:BlockSize:NumFinalTracks
% 		RightIDX = min(z + BlockSize - 1, NumFinalTracks);
% 		CurIDX = (z:RightIDX);
% 		FinalTracks(CurIDX) = tracks_world_to_img_trackvis(FinalTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
% 	end
% 	
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
	
% 	TracksHeader.init_threshold = 0.2000;
% 	TracksHeader.lmax = 6;
% 	TracksHeader.max_dist = 200;
% 	TracksHeader.max_num_attempts = NaN;
% 	TracksHeader.max_num_tracks = 5000;
% 	TracksHeader.max_trials = 50;
% 	TracksHeader.method = 'SD_PROB';
% 	TracksHeader.min_curv = 1;
% 	TracksHeader.min_dist = 10;
% 	TracksHeader.no_mask_interp = 0;
% 	TracksHeader.sh_precomputed = 1;
% 	TracksHeader.source = 'CSD.nii';
% 	TracksHeader.step_size = 0.2000;
% 	TracksHeader.stop_when_included = 0;
% 	TracksHeader.threshold = 0.1000;
% 	TracksHeader.unidirectional = 0;
% 	TracksHeader.roi.type = {'seed', 'exclude', 'mask'};
% 	TracksHeader.roi.file = {'freesurfer_wm.nii', 'freesurfer_wm_exclude.nii', 'mask.nii'};
% 
% 	TracksHeader.datatype.MATLABType = 'single';
% 	TracksHeader.datatype.Endianness = 'l';
% 	TracksHeader.count = 5000;
% 	TracksHeader.total_count = 22991;
% 
% 	TracksHeader.file = '. 439';
% 	TracksHeader.dataoffset = 439;
	
	%TracksDir = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
	%[~,~,~] = mkdir(TracksDir);
	%delete(fullfile(TracksDir, '*.tck'));
	%delete(fullfile(TracksDir, '*.trk'));
% 	for I = 1:length(CurLabels.shortlabels)
% 		for J = I + 1:length(CurLabels.shortlabels)
% 			T = ((FinalStartRegions == I) & (FinalEndRegions == J)) | ...
% 				((FinalStartRegions == J) & (FinalEndRegions == I));
% 			
% 			FIJ = ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_' CurLabels.shortlabels{I} '-' CurLabels.shortlabels{J}];
% % 			if(exist(fullfile(TracksDir, [FIJ '.tck']), 'file') == 2)
% % 				delete(fullfile(TracksDir, [FIJ '.tck']));
% % 				delete(fullfile(TracksDir, [FIJ '.trk']));
% % 			end
% % 			
% % 			FJI = ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_' CurLabels.shortlabels{J} '-' CurLabels.shortlabels{I}];
% % 			if(exist(fullfile(TracksDir, [FJI '.tck']), 'file') == 2)
% % 				delete(fullfile(TracksDir, [FJI '.tck']));
% % 				delete(fullfile(TracksDir, [FJI '.trk']));
% % 			end
% % 			
% 			if(any(T))
% 				TracksHeader.Tracks = FinalTracks(T);
% 				TracksHeader.count = length(TracksHeader.Tracks);
% 				
% 				save_mif_tracks(TracksHeader, fullfile(TracksDir, [FIJ '.tck']));% disp(['ln -sf -T ' FIJ '.tck' ' ' fullfile(TracksDir, [FJI '.tck'])]);
% 				[~, ~] = system(['ln -sf ' FIJ '.tck' ' ' fullfile(TracksDir, [FJI '.tck'])]);
% 				
% 				save_trackvis_tracks_mif(fullfile(TracksDir, [FIJ '.trk']), TracksHeader, AcceptedSeedsMIF);
% 				[~, ~] = system(['ln -sf ' FIJ '.trk' ' ' fullfile(TracksDir, [FJI '.trk'])]);
% 			end
% 			
% 		end
% 		F = fullfile(Subject, ...
% 			['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)], ...
% 			['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_' CurLabels.shortlabels{I}]);
% 		
% % 		if(exist([F '.tck'], 'file') == 2)
% % 			delete([F '.tck']);
% % 			delete([F '.trk']);
% % 		end
% % 		
% 		T = (FinalStartRegions == I) | (FinalEndRegions == I);
% 		if(any(T))	
% 			TracksHeader.Tracks = FinalTracks(T);
% 			TracksHeader.count = length(TracksHeader.Tracks);
% 			
% 			save_mif_tracks(TracksHeader, [F '.tck']);
% 			save_trackvis_tracks_mif([F '.trk'], TracksHeader, AcceptedSeedsMIF);	
% 		end
% 	end
	
% 	FinalRejectedTracks = tracks_img_to_world(FinalRejectedTracks, WMMaskMIF, size(WMMaskIMG));
% 	TracksHeader.Tracks = FinalRejectedTracks;
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_rejected_tracks_' num2str(Curvature)]);
%	save_mif_tracks(TracksHeader, [F '.tck']);
%	save_trackvis_tracks_mif([F '.trk'], TracksHeader, AcceptedSeedsMIF);
	%keyboard;
	SeedSizeMatrix = bsxfun(@plus, SeedSizes(:), SeedSizes(:)');
	SeedSizeMatrix = SeedSizeMatrix + double(SeedSizeMatrix == 0);
	SizeWeightedA = CountA ./ SeedSizeMatrix;
	WeightedA = WeightedA .* 2 ./ SeedSizeMatrix;
	T = CountA;
	T(T == 0) = 1;
	LengthA = LengthA ./ T;
	clear T;
	
	save(fullfile(Subject, ['connectivity_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']), 'FinalStartRegions', 'FinalEndRegions', 'WeightedA', 'CountA', 'SizeWeightedA', 'LengthA', 'SeedSizes');
	
	%FinalTracks = AcceptedTracks(TracksToAdd);
	
	%FinalAcceptedTracks = tracks_img_to_world(FinalAcceptedTracks, WMMaskMIF, size(WMMaskIMG));
	%FinalAcceptedTracks = AcceptedTracks;
% 	AcceptedTracks = tracks_img_to_world(AcceptedTracks, WMMaskMIF, size(WMMaskIMG));
% 	RejectedTracks = tracks_img_to_world(RejectedTracks, WMMaskMIF, size(WMMaskIMG));
% 	
 	AcceptedSeedsNII.img = int32(permute(flipdim(AcceptedSeedsIMG, 1), [2 1 3]));
 	save_nii(AcceptedSeedsNII, fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_seeds_' num2str(Curvature) '.nii.gz']));
% 	RejectedSeedsMIF.img = int32(permute(flipdim(RejectedSeedsIMG, 1), [2 1 3]));
% 	save_mif(RejectedSeedsMIF, fullfile(Subject, ['tracks_wb_eddy_cortex_' MrtrixMethod '_rejected_seeds_' num2str(Curvature) '.mif']));
% 	RejectedEndsMIF.img = int32(permute(flipdim(RejectedEndsIMG, 1), [2 1 3]));
% 	save_mif(RejectedSeedsMIF, fullfile(Subject, ['tracks_wb_eddy_cortex_' MrtrixMethod '_rejected_ends_' num2str(Curvature) '.mif']));
% 	
% 	% save culled tracks from anterior segment from CC
% 	TracksHeader.Tracks = AcceptedTracks(AcceptedAnteriorCCPosterior);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	F = fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_cc_anterior_posterior']);
% 	save_mif_tracks(TracksHeader, [F '.tck']);
% 	save_trackvis_tracks_mif([F '.trk'], TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = RejectedTracks(RejectedAnteriorCCPosterior);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	F = fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_cc_anterior_posterior']);
% 	save_mif_tracks(TracksHeader, [F '.tck']);
% 	save_trackvis_tracks_mif([F '.trk'], TracksHeader, AcceptedSeedsMIF);
% 	
% 	% save interhemispheric tracks that did not pass through the callosum
% % 	TracksHeader.Tracks = RejectedTracks(AnyInterNoCCRejected);
% % 	TracksHeader.count = length(TracksHeader.Tracks);
% % 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_no_cc.tck']));
% % 	TracksHeader.Tracks = AcceptedTracks(AnyInterNoCCAccepted);
% % 	TracksHeader.count = length(TracksHeader.Tracks);
% % 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_no_cc.tck']));
% 	
% 	% save all interhemispheric tracks
% % 	TracksHeader.Tracks = RejectedTracks(AnyInterRejected);
% % 	TracksHeader.count = length(TracksHeader.Tracks);
% % 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter.tck']));
% % 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter.trk']), TracksHeader, AcceptedSeedsMIF);
% % 	TracksHeader.Tracks = AcceptedTracks(AnyInterAccepted);
% % 	TracksHeader.count = length(TracksHeader.Tracks);
% % 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter.tck']));
% % 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = AcceptedTracks(AnyInterAnteriorCCAccepted);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_anterior_cc.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_anterior_cc.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = AcceptedTracks(~AcceptedAnteriorCCPosterior & AnyInterAnteriorCCAccepted);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_anterior_cc_culled.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_anterior_cc_culled.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = AcceptedTracks(AnyInterPosteriorCCAccepted);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_posterior_cc.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_inter_posterior_cc.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = RejectedTracks(~RejectedAnteriorCCPosterior & AnyInterAnteriorCCRejected);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_anterior_cc_culled.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_anterior_cc_culled.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = RejectedTracks(AnyInterAnteriorCCRejected);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_anterior_cc.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_anterior_cc.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = RejectedTracks(AnyInterPosteriorCCRejected);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_posterior_cc.tck']));
% 	save_trackvis_tracks_mif(fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '_inter_posterior_cc.trk']), TracksHeader, AcceptedSeedsMIF);
% 	
% 	TracksHeader.Tracks = RejectedTracks(~RejectedAnteriorCCPosterior);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_tracks_' num2str(Curvature) '.tck']));
% 	%FinalAcceptedTracks = tracks_img_to_world(FinalAcceptedTracks, WMMaskMIF, size(WMMaskIMG));
% 	TracksHeader.Tracks = AcceptedTracks(~AcceptedAnteriorCCPosterior);
% 	TracksHeader.count = length(TracksHeader.Tracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '.tck']));
% 	%AcceptedTracks = tracks_img_to_world(AcceptedTracks, WMMaskMIF, size(WMMaskIMG));
% 	%TracksHeader.Tracks = AcceptedTracks(InterhemisphericNoCC);
% 	%TracksHeader.count = length(TracksHeader.Tracks);
% 	%save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_interhemispheric_no_cc.tck']));
% 	
	%FinalTrackOneIMG = zeros(size(WMMaskIMG));
	
	
	%ISE = cellfun(@isempty, FinalTracks);
	%FinalTracks = FinalTracks(~ISE);
	%FinalTracks = tracks_img_to_world(FinalTracks, WMMaskMIF, size(WMMaskIMG));
% 	for z = 1:length(FinalTracks)
% 		%if(z == 1)
% % 			I = round(FinalTracks{z}(:, 2)) + size(FinalTrackOneIMG, 1) * round(FinalTracks{z}(:, 1)) + size(FinalTrackOneIMG, 1) * size(FinalTrackOneIMG, 2) * round(FinalTracks{z}(:, 3));
% % 			for k = 1:length(I)
% % 				FinalTrackOneIMG(I(k)) = FinalTrackOneIMG(I(k)) + 1;
% % 			end
% 		%end
% 		
% 		FinalTracks{z}(:, 1) = double(WMMaskMIF.dim(1)) - FinalTracks{z}(:, 1) + 1;
% 		FinalTracks{z}(:, 2) = double(WMMaskMIF.dim(2)) - FinalTracks{z}(:, 2) + 1;
% 		
% 		% change to 0-indexing
% 		FinalTracks{z} = FinalTracks{z} - 1;
% 		
% 		FinalTracks{z} = bsxfun(@times, FinalTracks{z}, WMMaskMIF.vox(1:3)');
% 		FinalTracks{z} = ScaledTransform * [FinalTracks{z}'; ones(1, size(FinalTracks{z}, 1))];
% 		FinalTracks{z} = FinalTracks{z}(1:3, :)';
% 	end
% 	FinalTrackOneMIF = WMMaskMIF;
% 	FinalTrackOneMIF.img = int32(permute(flipdim(FinalTrackOneIMG, 1), [2 1 3]));
% 	FinalTrackOneMIF.datatype.MATLABType = 'int32';
% 	save_mif(FinalTrackOneMIF, fullfile(Subject, ['tracks_' MrtrixMethod '_track_one_' num2str(Curvature) '.mif']));
% 	TracksHeader.Tracks = FinalTracks;
% 	TracksHeader.count = length(FinalTracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_' num2str(Curvature) '.tck']));
% 	
% 	TracksHeader.Tracks = FinalRejectedTracks;
% 	TracksHeader.count = length(FinalRejectedTracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_' num2str(Curvature) '.tck']));
 	
	%FinalTrackOneIMG = zeros(size(WMMaskIMG));
	
% 	for z = 1:length(FinalTracks)
% 		%if(z == 1)
% 			I = round(FinalTracks{z}(:, 2)) + size(FinalTrackOneIMG, 1) * round(FinalTracks{z}(:, 1)) + size(FinalTrackOneIMG, 1) * size(FinalTrackOneIMG, 2) * round(FinalTracks{z}(:, 3));
% 			
% 			for k = 1:length(I)
% 				FinalTrackOneIMG(I(k)) = FinalTrackOneIMG(I(k)) + 1;
% 			end
% 		%end
% 		
% 		FinalTracks{z}(:, 1) = double(WMMaskMIF.dim(1)) - FinalTracks{z}(:, 1) + 1;
% 		FinalTracks{z}(:, 2) = double(WMMaskMIF.dim(2)) - FinalTracks{z}(:, 2) + 1;
% 		
% 		% change to 0-indexing
% 		FinalTracks{z} = FinalTracks{z} - 1;
% 		
% 		FinalTracks{z} = FinalTracks{z} .* repmat(WMMaskMIF.vox(1:3)', size(FinalTracks{z}, 1), 1);
% 		FinalTracks{z} = ScaledTransform * [FinalTracks{z}'; ones(1, size(FinalTracks{z}, 1))];
% 		FinalTracks{z} = FinalTracks{z}(1:3, :)';
% 	end
% 	
	
% 	save(fullfile(Subject, ['tracks_' MrtrixMethod '_same_region']), 'SameRegionTracks');
% 	for z = 1:length(SameRegionTracks)
% 		SameRegionTracks{z}(:, 1) = double(WMMaskMIF.dim(1)) - SameRegionTracks{z}(:, 1) + 1;
% 		SameRegionTracks{z}(:, 2) = double(WMMaskMIF.dim(2)) - SameRegionTracks{z}(:, 2) + 1;
% 
% 		% change to 0-indexing
% 		SameRegionTracks{z} = SameRegionTracks{z} - 1;
% 
% 		SameRegionTracks{z} = ScaledTransform * [SameRegionTracks{z}'; ones(1, size(SameRegionTracks{z}, 1))];
% 		SameRegionTracks{z} = SameRegionTracks{z}(1:3, :)';
% 	end
% 
% 	TracksHeader.Tracks = SameRegionTracks;
% 	TracksHeader.count = length(SameRegionTracks);
% 	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_same_region.tck']));
% 	
%	RejectedTracksEndPointsMIF = WMMaskMIF;
%	RejectedTracksEndPointsMIF.datatype.MATLABType = 'int32';
%	RejectedTracksEndPointsIMG = zeros(size(WMMaskIMG), 'int32');
	%save(fullfile(Subject, ['tracks_' MrtrixMethod '_rejected']), 'RejectedTracks');
	
%	RejectedTracks = tracks_img_to_world(RejectedTracks, WMMaskMIF, size(WMMaskIMG));
%	RejectedTracksEndPointsIDX = zeros(1, length(RejectedTracks));
	
 	%for z = 1:length(RejectedTracks)
		%if(Codes(z) == -4 && TckSize(z) > MinSize(z))
	%		RejectedTracksEndPointsIDX(z) = round(RejectedTracks{z}(end, 2)) + size(RejectedTracksEndPointsIMG, 1) * round(RejectedTracks{z}(end, 1)) + size(RejectedTracksEndPointsIMG, 1) * size(RejectedTracksEndPointsIMG, 2) * round(RejectedTracks{z}(end, 3));
	%	end
		
 		%RejectedTracksEndPointsIMG(round(RejectedTracks{z}(1, 2)), round(RejectedTracks{z}(1, 1)), round(RejectedTracks{z}(1, 3))) = RejectedTracksEndPointsIMG(round(RejectedTracks{z}(1, 2)), round(RejectedTracks{z}(1, 1)), round(RejectedTracks{z}(1, 3))) + 1;
		%if(Codes(z) == -4)
		%	RejectedTracksEndPointsIMG(RejectedTracksEndPointsIDX(z)) = RejectedTracksEndPointsIMG(RejectedTracksEndPointsIDX(z);
		%end
 		
 	%	RejectedTracks{z}(:, 1) = double(WMMaskMIF.dim(1)) - RejectedTracks{z}(:, 1) + 1;
 	%	RejectedTracks{z}(:, 2) = double(WMMaskMIF.dim(2)) - RejectedTracks{z}(:, 2) + 1;
 
 		% change to 0-indexing
 	%	RejectedTracks{z} = RejectedTracks{z} - 1;
	%	RejectedTracks{z} = bsxfun(@times, RejectedTracks{z}, WMMaskMIF.vox(1:3)');
 	%	RejectedTracks{z} = ScaledTransform * [RejectedTracks{z}'; ones(1, size(RejectedTracks{z}, 1))];
 	%	RejectedTracks{z} = RejectedTracks{z}(1:3, :)';
	%end
	%N = histc(RejectedTracksEndPointsIDX(Codes == -4 & TckSize > MinSize), 1:numel(RejectedTracksEndPointsIMG));
	
	%RejectedTracksEndPointsIMG = reshape(N, size(RejectedTracksEndPointsIMG));
	
	%clear N;
	
	%[~, MaxI] = max(RejectedTracksEndPointsIMG(:));
	%[ID, JD, KD] = ind2sub(size(RejectedTracksEndPointsIMG), MaxI)
	%RejectedTracksEndPointsMIF.img = permute(flipdim(RejectedTracksEndPointsIMG, 1), [2 1 3]);
 	%save_mif(RejectedTracksEndPointsMIF, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_curvature_' num2str(Curvature) '.mif']));
	%RejectedTracksEndPointsIMG(:) = 0;
	%RejectedTracksEndPointsIMG(MaxI) = 1;
	%RejectedTracksEndPointsMIF.img = permute(flipdim(RejectedTracksEndPointsIMG, 1), [2 1 3]);
 	%save_mif(RejectedTracksEndPointsMIF, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_curvature_most_occurring_voxel_' num2str(Curvature) '.mif']));
	
 	%TracksHeader.Tracks = RejectedTracks(Codes == -4 & RejectedTracksEndPointsIDX == MaxI);
 	%TracksHeader.count = length(TracksHeader.Tracks);
 	%save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_most_occurring_voxel_' num2str(Curvature) '.tck']));

%  	TracksHeader.Tracks = RejectedTracks;
%  	TracksHeader.count = length(RejectedTracks);
%  	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_' num2str(Curvature) '.tck']));
% 	
% 	AcceptedThenRejectedTracks = tracks_img_to_world(AcceptedThenRejectedTracks, WMMaskMIF, size(WMMaskIMG));
% 	
% 	TracksHeader.Tracks = AcceptedThenRejectedTracks;
%  	TracksHeader.count = length(TracksHeader.Tracks);
%  	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_accepted_then_rejected_' num2str(Curvature) '.tck']));
% 	% save the different codes
	%UniqueCodes = unique(RejectedCodes);
	%for z = 1:length(UniqueCodes)
 	%	TracksHeader.Tracks = FinalRejectedTracks(RejectedCodes == UniqueCodes(z));
 	%	TracksHeader.count = length(TracksHeader.Tracks);
 	%	save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_code_' num2str(UniqueCodes(z)) '_' num2str(Curvature) '.tck']));
 	%end
	
% 	RejectedTracksEndPointsMIF.img = permute(flipdim(RejectedTracksEndPointsIMG, 1), [2 1 3]);
% 	save_mif(RejectedTracksEndPointsMIF, fullfile(Subject, ['tracks_' MrtrixMethod '_rejected_' num2str(Curvature) '.mif']));
	%[~,~,~] = mkdir(fullfile(Subject, ['tracks_' SeedType]));
% 	for CurLabel = 1:length(CurLabels.labels)
% 		M = (FinalStartRegions == CurLabel | FinalEndRegions == CurLabel);
% 		%keyboard;
% 		TracksHeader.Tracks = FinalTracks(M);
% 		TracksHeader.count = sum(M);
% 		save_mif_tracks(TracksHeader, fullfile(Subject, ['tracks' SeedType], ['tracks_' MrtrixMethod '_' num2str(CurLabels.values(CurLabel)) '.tck']));
% 	end
end

%%
% Degree = sum(CountA, 2);
% Degree = Degree ./ N;
% DegreeIMG = zeros(size(SeedIMG), 'int32');
% 
% for z = 1:length(CurLabels.values)
% 	DegreeIMG(SeedIMG == z) = int32(Degree(z));
% end
% 
% DegreeMIF = APARCSeedMIF;
% DegreeMIF.img = permute(flipdim(DegreeIMG, 1), [2 1 3]);
% DegreeMIF.datatype.MATLABType = 'int32';
% save_mif(DegreeMIF, fullfile(Subject, 'aparc_degree_img_' num2str(Curvature) '.mif'));
%%
% MIFTracks = load_mif_tracks(fullfile(Subject, 'tracks_wb_dt_stream.tck'));
% MIFTracks.count = sum(ValidTracks);
% MIFTracks.Tracks = MIFTracks.Tracks(ValidTracks);
% save_mif_tracks(MIFTracks, 'IRC143-1001/tracks_wb_dt_stream_culled.tck');
%return;
%%
% keyboard;
% I = find(StartRegions > 1000 & EndRegions > 1000);
% C = 13;
% subplot 121;
% hold off;
% imagesc(SeedIMG(:, :, round(CurTracks{I(C)}(1, 3))));
% hold on;
% plot(CurTracks{I(C)}(1, 1), CurTracks{I(C)}(1, 2), 'w*', 'MarkerSize', 20);
% subplot 122;
% hold off;
% imagesc(SeedIMG(:, :, round(CurTracks{I(C)}(end, 3))));
% hold on;
% plot(CurTracks{I(C)}(end, 1), CurTracks{I(C)}(end, 2), 'w*', 'MarkerSize', 20);
% %%
% return;
% for CurSeedTypeIDX = 3%1:length(SeedType)
% 	CurSeedType = SeedType{CurSeedTypeIDX};
% 	disp(CurSeedType);
% 	SeedIMG = zeros(size(WMMaskIMG), 'uint8');
% 	CountA = zeros(length(Labels{CurSeedTypeIDX}));
% 	WeightedA = zeros(length(Labels{CurSeedTypeIDX}));
% 	LengthA = zeros(length(Labels{CurSeedTypeIDX}));
% 	for CurLabel = 1:length(Labels{CurSeedTypeIDX})
% 		%fullfile(Subject, ['seeds_' CurSeedType], ['seed_' num2str(Labels{CurSeedTypeIDX}(CurLabel)) '_' num2str(Curvature) '.mif'])
% 		CurSeedMIF = load_mif(fullfile(Subject, ['seeds_' CurSeedType], ['seed_' num2str(Labels{CurSeedTypeIDX}(CurLabel)) '_' num2str(Curvature) '.mif']));
% 		CurSeedIMG = flipdim(permute(CurSeedMIF.img, [2 1 3:ndims(CurSeedMIF.img)]), 1);
% 		SeedIMG(CurSeedIMG > 0) = CurLabel;
% 		clear CurSeedIMG CurSeedMIF;
% 	end
% 	SeedSizes = histc(SeedIMG(SeedIMG > 0), 1:length(Labels{CurSeedTypeIDX}));
% 	save(fullfile(Subject, ['seed_sizes_' CurSeedType '.mat']), 'SeedSizes');
% 	%clear SeedSizes;
% 	SeedMIF = WMMaskMIF;
% 	SeedMIF.img = permute(flipdim(SeedIMG, 1), [2 1 3]);
% 	save_mif(SeedMIF, fullfile(Subject, ['seeds_' CurSeedType '.mif']));
% 	%keyboard;
% AcceptedThenRejectedTracks
% 	TrackArcLengths = cellfun(@arc_length, CurTracks);
% 
% 	TracksSZ = cellfun(@(x) (size(x, 1)), CurTracks);
% 
% 	CurTracks = cat(1, CurTracks{:});
% 	CurTracks = InvTransform * [CurTracks'; ones(1, size(CurTracks, 1))];
% 	CurTracks = CurTracks(1:3, :)';
% 
% 	% change to 1-indexing
% 	CurTracks = CurTracks + 1;
% 
% 	CurTracks(:, 1) = double(WMMaskMIF.dim(1)) - CurTracks(:, 1) + 1;
% 	CurTracks(:, 2) = double(WMMaskMIF.dim(2)) - CurTracks(:, 2) + 1;
% 
% 	CurTracks = mat2cell(CurTracks, TracksSZ, 3);
% 	CurTrackEndRegions = zeros(1, length(CurTracks), 'uint8');
% 	WeightedConnectivityCounts = zeros(1, length(Labels{CurSeedTypeIDX}));
% 	LengthConnectivityCounts = zeros(1, length(Labels{CurSeedTypeIDX}));
% 	for CurTrackIDX = 1:length(CurTracks)
% 		ID = int32(round(CurTracks{CurTrackIDX}(:, 2)));
% 		JD = int32(round(CurTracks{CurTrackIDX}(:, 1)));
% 		KD = int32(round(CurTracks{CurTrackIDX}(:, 3)));
% 		%I = sub2ind(size(SeedIMG), ID, JD, KD);
% 		% inline sub2ind
% 		I = ID + int32(WMMaskIMGSZ(1)) * (JD - 1) + int32(NumInSlice) * (KD - 1);
% 		CurRegions = SeedIMG(I);
% 		CurRegions = CurRegions(CurRegions > 0);
% 
% 		if(~isempty(CurRegions))
% 			% if the last value is the current seed then swap the 
% 			if(CurRegions(end) == CurLabel)
% 				if(CurRegions(1) ~= CurLabel)
% 					CurTrackEndRegions(CurTrackIDX) = CurRegions(1);
% 				end
% 			else
% 				CurTrackEndRegions(CurTrackIDX) = CurRegions(end);
% 			end
% 			if(CurTrackEndRegions(CurTrackIDX) > 0)
% 				WeightedConnectivityCounts(CurTrackEndRegions(CurTrackIDX)) = WeightedConnectivityCounts(CurTrackEndRegions(CurTrackIDX)) + 1 ./ TrackArcLengths(CurTrackIDX);
% 				LengthConnectivityCounts(CurTrackEndRegions(CurTrackIDX)) = LengthConnectivityCounts(CurTrackEndRegions(CurTrackIDX)) + TrackArcLengths(CurTrackIDX);
% 				TrackLengthsByROI{CurLabel} = [TrackLengthsByROI{CurLabel} TrackArcLengths(CurTrackIDX)];
% 				TrackLengthsByROI{CurTrackEndRegions(CurTrackIDX)} = [TrackLengthsByROI{CurTrackEndRegions(CurTrackIDX)} TrackArcLengths(CurTrackIDX)];
% 			end
% 		end
% 	end
% 	clear ID JD KD I CurRegions;
% 	%CurTracks = cellfun(@(x) (int16(round(x))), CurTracks, 'UniformOutput', false);
% 	ConnectivityCounts = hiAcceptedThenRejectedTracksstc(double(CurTrackEndRegions), 1:length(Labels{CurSeedTypeIDX}));
% 
% 	CountA(CurLabel, :) = ConnectivityCounts;
% 	WeightedA(CurLabel, :) = WeightedConnectivityCounts;
% 	LengthA(CurLabel, :) = LengthConnectivityCounts;
% end
% 
% SeedSizeMatrix = bsxfun(@plus, SeedSizes, SeedSizes');
% CountA = CountA + CountA';
% WeightedA = 2 .* (WeightedA + WeightedA') ./ SeedSizeMatrix;
% B = CountA;
% B(CountA == 0) = 1;
% LengthA = (LengthA + LengthA') ./ B;
% LengthA(CountA == 0) = Inf;
% clear B;
% 
% %	save(fullfile(Subject, ['connectivity_' CurSeedType '.mat']), 'CountA', 'WeightedA', 'LengthA', 'TrackLengthsByROI');
% %end
% 
% % 	CurTrackEndPoints = cellfun(@(x) (x(end, :)), CurTracks, 'UniformOutput', false);
% % 	CurTrackEndPoints = cell2mat(CurTrackEndPoints);
% % 
% % 	ID = round(CurTrackEndPoints(:, 2));
% % 	JDChallenges and limitations of quantifying brain connectivity in vivo with diffusion MRI = round(CurTrackEndPoints(:, 1));
% % 	KD = round(CurTrackEndPoints(:, 3));
% % 	I = sub2ind(size(IncludeIMG), ID, JD, KD);
% % 	
% % 	CurEndRegions = SeedIMG(I);
% % 	F = find(CurEndRegions == 1, 1, 'first');
% % 	clf;
% % 	imshow(IncludeIMG(:, :, KD(F)), []);
% % 	hold on;
% % 	plot(JD(F), ID(F), 'r*');
% 	%N = histc(CurEndRegions, [0:length(CurLabels)]);
% 	%keyboard;
% 	%end
% % 	ID = round(CurTracks(:, 2));
% % 	JD = round(CurTracks(:, 1));
% % 	KD = round(CurTracks(:, 3));
% 	
% 	
% 	%%
% 	
% 	%%
% 
% 
% 
% % FreesurferLabels = read_freesurfer_ctx_labels(fullfile('..', 'FreeSurferColorLUT.txt'));
% % 
% % FreesurferValues = [FreesurferLabels.value];
% % [TF, LOC] = ismember(CurLabels, FreesurferValues);
% % [FL{1:length(FreesurferLabels)}] = deal(FreesurferLabels.label);
% % D = FL(LOC);
% % %%
% % axis on;
% % set(gca, 'XTick', 1:68, 'YTick', 1:68, 'XTickLabel', [], 'YTickLabel', [])
% % 
% % text(1:length(A), repmat(-1, 1, length(A)), D, 'rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% % text(repmat(-1, 1, length(A)), 1:length(A), D, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
% 
% % 
% % APARCNII = load_nii(fullfile(Subject, 'aparc+aseg.nii.gz'));
% % APARCIMG = flipdim(permute(APARCNII.img, [2 1 3]), 1);
% % APARCNII.img = [];
% % APARCBANII = load_nii(fullfile(Subject, 'aparc.a2009s+aseg.nii.gz'));
% % APARCBAIMG = flipdim(permute(APARCBANII.img, [2 1 3]), 1);
% % APARCBANII.img = [];
% % 
% % %aparc+aseg 1000's are lh cortex
% % %aparc+aseg 2000's are rh cortex
% % 
% % %aparc.2009a+aseg 11000's are lh cortex
% % %aparc.2009a+aseg 12000's are rh cortex
% % 
% % FreesurferLabels = read_freesurfer_ctx_labels('FreeSurferColorLUT.txt');
% % FreesurferValues = cat(1, FreesurferLabels.value);
% % APARCValues = FreesurferValues(FreesurferValues >= 1000 & FreesurferValues < 3000);
% % load(fullfile(Subject, 'whole_brain_tracts'));
% % % convert to single to save memory
% % 
% % RegVoxTracts = cellfun(@single, RegVoxTracts, 'UniformOutput', false);
% % TractSZ = cellfun(@(x) (size(x, 1)), RegVoxTracts);
% % 
% % RegVoxTracts = cat(1, RegVoxTracts{:});
% % 
% % I = RegVoxTracts(:, 2);
% % J = RegVoxTracts(:, 1);
% % K = RegVoxTracts(:, 3);
% % 
% % ID = round(I); JD = round(J); KD = round(K);
% % clear I J K;
% % IDX = sub2ind(size(APARCIMG), ID, JD, KD);
% % clear ID JD KD;
% % LabelsSeen = APARCIMG(IDX);
% % clear IDX;
% % LabelsSeen = mat2cell(LabelsSeen, TractSZ, 1);
% % LabelsSeen = cellfun(@unique, LabelsSeen, 'UniformOutput', false);
% % LabelsSeen = cellfun(@(x) (intersect(x, APARCValues)), LabelsSeen, 'UniformOutput', false);
% % 
% % %LabelsSeenIDX = LabelsSeen;
% % %%
% % A = zeros(length(APARCValues), length(APARCValues));
% % 
% % for z = 1:length(LabelsSeen)
% % 	[TF, LOC] = ismember(LabelsSeen{z}, APARCValues);
% % 	%LabelsSeenIDX{z} = LOC;
% % 	if(length(LOC) > 2)
% % 		[IIDX, JIDX] = meshgrid(LOC, LOC);
% % 		I = sub2ind(size(A), IIDX(:), JIDX(:));
% % 		A(I) = A(I) + 1;
% % 	end
% % end
% % A = A .* (1 - eye(length(A)));
% % 
% % % ID = floor(I); JD = floor(J); KD = floor(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = ceil(I); JD = floor(J); KD = floor(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = floor(I); JD = ceil(J); KD = floor(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = floor(I); JD = floor(J); KD = ceil(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = ceil(I); JD = ceil(J); KD = floor(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = ceil(I); JD = floor(J); KD = ceil(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = floor(I); JD = ceil(J); KD = ceil(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % % ID = ceil(I); JD = ceil(J); KD = ceil(K); IDX = sub2ind(size(BZeroSwappedIMG), ID, JD, KD);
% % 
% % %%
% % %subplot 121;
% % % clf;
% % % H = slice(double(APARCIMG), size(APARCIMG, 2) / 2, size(APARCIMG, 1) / 2, size(APARCIMG, 3) / 2 + -6);
% % % set(H, 'EdgeColor', 'none', 'AmbientStrength', 1);
% % % colormap gray;
% % % hold on;
% % % streamline(RegVoxTracts(1:1000));
% % % axis ij equal;
% % % view(2);
% % % xlabel('x');
% % % ylabel('y');
% % % zlabel('z');
% % 

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

% for k = 1:length(SeedValues)
% 	if(any(SeedValues{k} == 255))
% 		% we will cull if any of the start or end points are posterior of
% 		% the central CC centroid
% 		if(any(S{k}(:, 2) > CentralCCCentroid))
% 			M(k) = 1;
% 		end
% 	end
% end
%clear SeedValues IDX TracksSZ;


% T = cat(1, S{:});
% TracksSZ = cellfun('size', S, 1);
% 
% SeedValues = interp3_linear_fast(double(SeedIMG), double(T(:, 1)), double(T(:, 2)), double(T(:, 3)), 'nearest');
% SeedValues = mat2cell_vec(SeedValues, TracksSZ);
% clear T;

% for z = 1:numel(S)
% 	% 255 is the label for the anterior segment
% 	if(any(SeedValues{z} == 255))
% 		% we will cull if any of the start or end points are posterior of
% 		% the central CC centroid
% 		if(any(S{z}(:, 2) > CentralCCCentroid))
% 			M(z) = 1;
% 		end
% 	end
% end
	
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
