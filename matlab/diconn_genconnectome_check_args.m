function [CurLabels, IncludeFile, ExcludeFile] = diconn_genconnectome_check_args(Subject, MrtrixMethod, SeedType, Curvature)

% mrtrix_create_tracks_check_args(Subject, MrtrixMethod, SeedType, Curvature)
% DESCRIPTION
% 	Helper function for mrtrix_create_tracks* functions to check the arguments
%	to see if we are choosing a valid MrtrixMethod, SeedType, Curvature and if all the files are there.
%	We also return the CurLabels, IncludeFile, ExcludeFile.
%
%ConnectomeSeedDir = mrtrix_create_tracks_mrtrixdir;
% MrtrixDir = 'mrtrix2DWISpaceNoBReorient';
ConnectomeSeedDir = 'ConnectomeSeedImages';

CortexLabels = load_freesurfer_cortex_labels('/usr/local/diconnectivity/lib/parc_schemes');
%CortexLabels = load_freesurfer_cortex_labels('parc_schemes');

L =  {'dt_stream', 'sd_stream', 'sd_prob'};
%keyboard;
if(~ismember(MrtrixMethod, L))
	error(['MrtrixMethod must be one of: ' sprintf('%s ', L{1:end-1}), L{end}]);
end
%keyboard;
switch(lower(SeedType))
	case 'aparc.a2009s'
		CurLabels = CortexLabels.APARCOhNine;
		ExcludeFile = fullfile(ConnectomeSeedDir, Subject, 'excludes_aparcohnine');
		IncludeFile = fullfile(ConnectomeSeedDir, Subject, 'seeds_aparcohnine');
	case 'aparc.a2005s'
		CurLabels = CortexLabels.APARCOhFive;
		ExcludeFile = fullfile(ConnectomeSeedDir, Subject, 'excludes_aparcohfive');
		IncludeFile = fullfile(ConnectomeSeedDir, Subject, 'seeds_aparcohfive');
	case 'subdivided1000'
		CurLabels = CortexLabels.SUBDIVIDED;
		ExcludeFile = fullfile(ConnectomeSeedDir, Subject, ['excludes_' SeedType ]);
		IncludeFile = fullfile(ConnectomeSeedDir, Subject, ['seeds_' SeedType ]);
	case 'wmsubdivided1000'
		CurLabels = CortexLabels.WMSUBDIVIDED;
		ExcludeFile = fullfile(ConnectomeSeedDir, Subject, ['excludes_' SeedType ]);
		IncludeFile = fullfile(ConnectomeSeedDir, Subject, ['wm_' SeedType ]);
	otherwise
		if(isfield(CortexLabels, upper(SeedType)))
			CurLabels = CortexLabels.(upper(SeedType));
			ExcludeFile = fullfile(ConnectomeSeedDir, Subject, ['excludes_' lower(SeedType) ]);
			IncludeFile = fullfile(ConnectomeSeedDir, Subject, ['seeds_' lower(SeedType) ]);
		else
			error(['Unsupported seed scheme: ' SeedType]);
		end
end

% FilesToCheck = {...
% 	IncludeFile, ...
% 	ExcludeFile, ...
% 	fullfile(MrtrixDir, Subject, 'dwi'), ...
% 	fullfile(ConnectomeSeedDir, Subject, 'freesurfer_wm'), ...
% 	fullfile(MrtrixDir, Subject, 'mask'), ...
% 	fullfile(MrtrixDir, Subject, 'grad.b'), ...
% 	fullfile(MrtrixDir, Subject, 'CSD'), ...
% };
% 
% for z = 1:length(FilesToCheck)
% 	if(exist(FilesToCheck{z}, 'file') ~= 2 && exist([FilesToCheck{z} '.nii.gz'], 'file') ~= 2)
% 		error(['File does not exist, rerun any required preprocessing steps: ' FilesToCheck{z}]);
% 	end
% end

if(isscalar(Curvature) && isnumeric(Curvature))
	if(Curvature < 0)
		error('Curvature must be positive');
	end
else
	error('Curvature must be a numeric scalar');
end
