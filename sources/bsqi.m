function [ N_J, N_G ] = bsqi( ann1, ann2, THR )

if nargin < 3
    THR = 0.15; % window for matching two peaks
end

% Use histc to bin each jqrs into bins centered on gqrs
xi = [ann1' - THR;
    ann1' + THR];
xi = xi(:);

% we now replace monotonically non-decreasing entries with their avg
% normally, we have two ECG beats next to each other, and windows do
% not overlap, i.e.:   |  x  |     l  x  l

% where 'x' is a beat with || and ll define the +- window

% sometimes the beats are too close however, and the windows overlap
% the following code changes windows like this:  |  l |  l
% to windows like this:                          |   |l  l

idxFix = [false;diff(xi) < 0];
xi_fixed = [xi(idxFix),xi([idxFix(2:end);false])];
xi_fixed = mean(xi_fixed,2);

xi(idxFix) = xi_fixed;
xi([idxFix(2:end);false]) = xi_fixed;

%=== now place gqrs into bins centered on jqrs
N_J = histc(ann2,xi);

%=== repeat for the other way round
xi = [ann2' - THR;
    ann2' + THR];
xi = xi(:);
idxFix = [false;diff(xi) < 0];
xi_fixed = [xi(idxFix),xi([idxFix(2:end);false])];
xi_fixed = mean(xi_fixed,2);

xi(idxFix) = xi_fixed;
xi([idxFix(2:end);false]) = xi_fixed;
N_G = histc(ann1,xi);


% for each peak in gqrs, N_G indicates whether it has a matching peak
% in jqrs (and vice versa for N_J)
N_J = N_J(1:2:end);
N_G = N_G(1:2:end);
N_J = N_J(:);
N_G = N_G(:);

end