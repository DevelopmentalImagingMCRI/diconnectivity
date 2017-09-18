function [FornixTracks, AnyInterhemispheric, AnyInterhemisphericNoCC, AnyInterhemisphericAnteriorCC, AnyInterhemisphericPosteriorCC] = mrtrix_create_tracks_through_fornix(S, SeedIMG)

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
