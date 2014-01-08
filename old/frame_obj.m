function [fx,grad] = frame_obj(data,regmax,pos,shape,template,eta,sig,xi)
% Gives the log probability of the data and regularized maxima for a given
% position and unnormalized shape (or shape*intesity)
%
% eta - observation noise
% sig - position noise
% xi  - difference-from-template "noise"

assert(~any(size(regmax)~=size(pos)));
assert(~any(size(data)~=size(template)));

