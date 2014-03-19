function J = cnormalize( I, sigSub, sigDiv  )
% Contrast-Normalize an image
%
% Input:
%   I - the image to normalize
%   sigSub - the kernel width for subtractive normalization
%   sigDiv - the kernel width for divisive normalization
%
% David Pfau, 2014

J = I - blur( I, sigSub );
J = J ./ blur( J.^2 - blur( J, sigDiv ).^2, sigDiv );