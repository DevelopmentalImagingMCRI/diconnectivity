function [CurLabels, IncludeFile, ExcludeFile] = mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature)

% mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature)
% DESCRIPTION
% 	Helper function for mrtrix_create_tracks* functions to check the arguments
%	to see if we are choosing a valid MrtrixMethod, SeedType, Curvature and if all the files are there.
%	We also return the CurLabels, IncludeFile, ExcludeFile.
%
MrtrixSubjDir = mrtrix_create_tracks_mrtrixdir;
MrtrixSubjDir = 'mrtrix2DWISpaceBReorient';
CortexLabels = load_freesurfer_cortex_labels;
L =  {'dt_stream', 'sd_stream', 'sd_prob'};
%keyboard;
if(~ismember(MrtrixMethod, L))
	error(['MrtrixMethod must be one of: ' sprintf('%s ', L{1:end-1}), L{end}]);
end
%keyboard;
switch(lower(SeedType))
	case 'aparc.a2009s'
		CurLabels = CortexLabels.APARCOhNine;
		ExcludeFile = fullfile(MrtrixSubjDir, Subject, 'exclude_aparcohnine.nii.gz');
		IncludeFile = fullfile(MrtrixSubjDir, Subject, 'brain_non_wm_aparcohnine.nii.gz');
	case 'aparc.a2005s'
		CurLabels = CortexLabels.APARCOhFive;
		ExcludeFile = fullfile(MrtrixSubjDir, Subject, 'exclude_aparcohfive.nii.gz');
		IncludeFile = fullfile(MrtrixSubjDir, Subject, 'brain_non_wm_aparcohfive.nii.gz');
	case 'subdivided1000'
		CurLabels = CortexLabels.SUBDIVIDED;
		ExcludeFile = fullfile(MrtrixSubjDir, Subject, ['exclude_' SeedType '.nii.gz']);
		IncludeFile = fullfile(MrtrixSubjDir, Subject, ['brain_non_wm_' SeedType '.nii.gz']);
	case 'wmsubdivided1000'
		CurLabels = CortexLabels.WMSUBDIVIDED;
		ExcludeFile = fullfile(MrtrixSubjDir, Subject, ['exclude_' SeedType '.nii.gz']);
		IncludeFile = fullfile(MrtrixSubjDir, Subject, ['wm_' SeedType '.nii.gz']);
	otherwise
		if(isfield(CortexLabels, upper(SeedType)))
			CurLabels = CortexLabels.(upper(SeedType));
			ExcludeFile = fullfile(MrtrixSubjDir, Subject, ['exclude_' lower(SeedType) '.nii.gz']);
			IncludeFile = fullfile(MrtrixSubjDir, Subject, ['brain_non_wm_' lower(SeedType) '.nii.gz']);
		else
			error(['Unsupported seed scheme: ' SeedType]);
		end
end

FilesToCheck = {...
	IncludeFile, ...
	ExcludeFile, ...
	fullfile(MrtrixSubjDir, Subject, 'dwi.nii.gz'), ...
	fullfile(MrtrixSubjDir, Subject, 'freesurfer_wm.nii.gz'), ...
	fullfile(MrtrixSubjDir, Subject, 'mask.nii.gz'), ...
	fullfile(MrtrixSubjDir, Subject, 'grad.b'), ...
	fullfile(MrtrixSubjDir, Subject, 'CSD.nii.gz'), ...
};

for z = 1:length(FilesToCheck)
	if(exist(FilesToCheck{z}, 'file') ~= 2)
		error(['File does not exist, rerun any required preprocessing steps: ' FilesToCheck{z}]);
	end
end

if(isscalar(Curvature) && isnumeric(Curvature))
	if(Curvature < 0)
		error('Curvature must be positive');
	end
else
	error('Curvature must be a numeric scalar');
end
