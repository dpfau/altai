function x = sample_logliks(logliks)

liks = exp(logliks-max(logliks)); % Prevents numerical overflow when converting from log domain. Underflow isn't a problem though.
cdf = cumsum(liks)/sum(liks);
x = find(cdf>rand,1);