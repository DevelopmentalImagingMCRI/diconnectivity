function [CortexLabels] = load_freesurfer_cortex_labels(varargin)

if nargin == 1
	InputDir = varargin{1};
else
	InputDir = '/home/addo/connectivity/trunk/dfbiconnectivity/lib/parc_schemes';
	InputDir = '/usr/local/diconnectivity/lib/parc_schemes';
end

if exist(InputDir, 'dir') ~= 7
	error(['Input directory does not exist: ' InputDir]);
end
D = dir(InputDir);

SchemeNameRemapSource{1} = 'APARC_ADS';
SchemeNameRemapDest{1} = 'APARCADS';

for z = 1:length(D)
	%disp([D(z).name ' ' num2str(exist(fullfile(InputDir, D(z).name), 'file'))]);
	if(exist(fullfile(InputDir, D(z).name), 'file') == 2)
		SchemeName = D(z).name;
		StructName = SchemeName;
% 		[TF, LOC] = ismember(SchemeName, SchemeNameRemapSource);
% 		if(TF)
% 			StructName = SchemeNameRemapDest{LOC};
% 		else
% 			StructName = SchemeName;
% 		end
		
		% skip schemes that have invalid field names
		% A valid MATLAB variable name is a character string of letters,
		% digits, and underscores, such that the first character is a letter,
		% and the length of the string is less than or equal to the value returned
		% by the namelengthmax function. Any string that exceeds namelengthmax is
		% truncated in the varname output. See Example 6, below.
		
		[match] = regexp(SchemeName, '[a-zA-Z_][a-zA-Z0-9_]*', 'match');
		if(isempty(match))
			warning(['File ' SchemeName ' has an invalid name for MATLAB']);
		else
			FID = fopen(fullfile(InputDir, D(z).name), 'r');
			tline = fgetl(FID);
			if(~ischar(tline))
				warning(['File ' D(z).name ' in the parc_schemes directory is empty']);
			else
				[match, tokens] = regexp(tline, '^Seed: (\S+)$', 'match', 'tokens');
				if(isempty(match))
					warning(['File ' D(z).name ' in the parc_schemes directory has an invalid first line']);
				else
					CortexLabels.(StructName).SeedFile = tokens{1}{1};
					A = textscan(FID, '%s %s %s');
					
					CortexLabels.(StructName).values = cellfun(@comma_split, A{1}, 'UniformOutput', false);
					%keyboard;
					%CortexLabels.(StructName).values = textscan(A{1}, '%d', 'Delimiter', ',');
					%keyboard;
					%CortexLabels.(StructName).values = CortexLabels.(StructName).values{1};
					%$#for k = 1:length(CortexLabels.(StructName).values)
					%	CortexLabels.(StructName).values{z} = textscan(CortexLabels.(StructName).values{z}, '%d,');
					%end
					CortexLabels.(StructName).labels = A{2};
					CortexLabels.(StructName).shortlabels = A{3};
					clear A;
				end
			end
			fclose(FID);
		end
	end
end

function [C] = comma_split(A)

C = textscan(A, '%d', 'Delimiter', ',', 'CollectOutput', true);
%C = C{1};
