function [ d, beat, q, p ] = resampleSignal(d, fsorig, fsdown, beat)
%RESAMPLESIGNAL	Downsamples signal d from fs2 Hz to fs1 Hz
%	[ d ] = resampleSignal(d, fsorig, fsdown) uses rat to find integer
%	approximations to resampling signal d at fsdown.
%
%
%	Inputs:
%       d   - vector data
%       fsorig - original sampling frequency
%       fsdown - downsampled frequency
%
%	Outputs:
%       d - downsampled signal
%
%
%	Example
%		[ ] = resampleSignal()
%
%	See also FCN1

%	References:
%
%

%	Copyright 2012 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 20-Jul-2012 14:59:53
%	Contact: alistairewj@gmail.com
[p,q] = rat(fsdown/fsorig,0.0001);
d = resample(d,p,q); % Integer resample with alias filter

if nargin>3
    beat = round(beat / q * p);
else
    beat = [];
end
if fsdown>fsorig
    fprintf('Warning: upsampling, artefacts possible.\n');
end

end