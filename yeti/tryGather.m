function x = tryGather(x)
% An idiom we use a lot to make the code work whether x is a gpuArray or
% not

try gpuDevice;
    x = gather(x);
catch e
end