function M = tracks_interhemispheric_no_cc(S, SeedIMG)

% marks the interhemispheric tracks that do not go through the callosum or
% through the white matter as 1

M = false(size(S));

T = cat(1, S{:});
TracksSZ = cellfun(@(x) (size(x, 1)), S);

SeedValues = interp3(double(SeedIMG), double(T(:, 1)), double(T(:, 2)), double(T(:, 3)), 'nearest');
SeedValues = mat2cell_vec(SeedValues, TracksSZ);
clear T;

for z = 1:length(S)
	LWM = (SeedValues{z} == 2);
	RWM = (SeedValues{z} == 41);
	
	CCLabels = [251 252 253 254 255];
	if(~any(...
			(RWM(2:end) & LWM(1:end - 1)) | ...
			(LWM(2:end) & RWM(1:end - 1)) ...
			) && all(~ismember(SeedValues{z}, CCLabels)))
		%keyboard;
		M(z) = 1;
	end
end