function A = arc_length(V)

% A = arc_length(V)
%
% DESCRIPTION
%	Computes the arc length A of the curve given by the vertices V. V can
%	have any arbitrary number of dimensions.
%	Arc length is computed as:
%	sum(i = 1 ... N - 1)||V(i) - V(i + 1)||_2
%
% PARAMETERS
%	V [NumVertices, NumDims]: the ordered vertices that form a contour
%		in Euclidean NumDims-space
%
% RETURNS
%	A [1]: the arc length

if(size(V, 1) > 1)
	VDash = diff(V, 1, 1);

	MagVDash = sqrt((sum(VDash .* VDash, 2)));

	A = double(sum(MagVDash, 1));
else
	A = double(0);
end