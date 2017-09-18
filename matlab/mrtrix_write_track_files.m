function [varargout] = mrtrix_write_track_files(Subject, MrtrixMethod, SeedType, Curvature, WhichTracks, TrackFileFormat)

% mrtrix_create_tracks_wb_wm(Subject, MrtrixMethod, SeedType, Curvature, WhichTracks, TrackFileFormat)
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
%	WhichTracks (char):
%		'accepted': writes all accepted tracks in the subject's directory
%		'rois': writes tracks for each roi 
% OUTPUT FILES

[CurLabels, IncludeFile, ExcludeFile] = mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature);

TrackFile = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_curvature_' num2str(Curvature) '.mat']);

if(~ismember(lower(WhichTracks), {'accepted', 'rois'}))
	error('WhichTracks is not supported');
end

if(~ismember(lower(TrackFileFormat), {'mrtrix', 'trackvis'}))
	error('TrackFileFormat not supported');
end
		
if(exist(TrackFile, 'file') ~= 2)
	error(['Track file: ' TrackFile ' is not a regular file']);
else

	load(TrackFile, 'FinalTracks', 'FinalArcLengths', 'FinalStartRegions', 'FinalEndRegions');
	NumFinalTracks = numel(FinalTracks);
	WMMaskNII = load_nii(fullfile(Subject, 'freesurfer_wm'));
	WMMaskIMG = flipdim(permute(WMMaskNII.img, [2 1 3]), 1);

	OriginalSeedNII = load_nii(fullfile('..', 'RegisterFreesurferToMrtrix', Subject, CurLabels.SeedFile));
	OriginalSeedIMG = flipdim(permute(OriginalSeedNII.img, [2 1 3]), 1);
	SeedNII = load_nii(IncludeFile);
	SeedIMG = flipdim(permute(SeedNII.img, [2 1 3]), 1);

	% tracks are in IMAGE coordinates

	EmptyTracksHeader.init_threshold = 0.2000;
	EmptyTracksHeader.lmax = 8;
	EmptyTracksHeader.max_dist = 200;
	EmptyTracksHeader.max_num_attempts = 200000000;
	EmptyTracksHeader.max_num_tracks = 2000000;
	EmptyTracksHeader.max_trials = 50;
	EmptyTracksHeader.method = upper(MrtrixMethod);
	EmptyTracksHeader.min_curv = Curvature;
	EmptyTracksHeader.min_dist = 10;
	EmptyTracksHeader.no_mask_interp = 0;
	EmptyTracksHeader.sh_precomputed = 1;
	EmptyTracksHeader.source = '/tmp/tmp.ndDE2n.nii.gz';
	EmptyTracksHeader.step_size = 1;
	EmptyTracksHeader.stop_when_included = 1;
	EmptyTracksHeader.threshold = 0.1000;
	EmptyTracksHeader.unidirectional = 0;
	EmptyTracksHeader.roi.type = {'seed', 'include', 'exclude', 'mask'};
	EmptyTracksHeader.roi.file = {...
		fullfile(Subject, 'freesurfer_wm.nii.gz'), ...
		IncludeFile, ...
		ExcludeFile, ...
		fullfile(Subject, 'mask.nii.gz')};
	EmptyTracksHeader.datatype.MATLABType = 'single';
	EmptyTracksHeader.datatype.Endianness = 'l';
	EmptyTracksHeader.count = NumFinalTracks;
	EmptyTracksHeader.total_count = NumFinalTracks;
	EmptyTracksHeader.file = '. 637';
	EmptyTracksHeader.dataoffset = 637;
	%keyboard;
	% no longer are we saving track files automatically, do it on request
	
	%% block to save RAM
	
	%T = tracks_img_to_world(FinalTracks, WMMaskNII, size(WMMaskIMG));
	%clear T;
	BlockSize = 50000;
	
	for z = 1:BlockSize:NumFinalTracks
		RightIDX = min(z + BlockSize - 1, NumFinalTracks);
		CurIDX = (z:RightIDX);
	
		FinalTracks(CurIDX) = tracks_img_to_world(FinalTracks(CurIDX), WMMaskNII, size(WMMaskIMG));
	end
	
	TracksHeader = EmptyTracksHeader;
	
	switch(lower(WhichTracks))
		case 'accepted'
			TracksHeader.Tracks = FinalTracks;
			TracksHeader.count = length(TracksHeader.Tracks);
			F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature)]);
			switch(lower(TrackFileFormat))
				case 'mrtrix'
					save_mif_tracks(TracksHeader, [F '.tck']);
				case 'trackvis'
					save_trackvis_tracks_nii([F '.trk'], TracksHeader, WMMaskNII);
			end
		case 'rois'
			[~,~,~] = mkdir(fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois']));
			for I = 1:length(CurLabels.shortlabels)
				TracksHeader.Tracks = FinalTracks((FinalStartRegions == I) | (FinalEndRegions == I));
				F = fullfile(Subject, ['tracks_' SeedType '_' MrtrixMethod '_accepted_tracks_' num2str(Curvature) '_rois'], CurLabels.shortlabels{I});
				switch(lower(TrackFileFormat))
					case 'mrtrix'
						save_mif_tracks(TracksHeader, [F '.tck']);
					case 'trackvis'
						save_trackvis_tracks_nii([F '.trk'], TracksHeader, WMMaskNII);
				end
			end
	end
end
