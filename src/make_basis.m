function basis = make_basis(mode,order,window)

basis = zeros(window,window,0);
switch mode
    case 'hermite-gaussian'
        x = linspace(-5,5,window);
        hg = @(n) hermite(n,x).*exp(-x.^2/2)/norm(hermite(n,x).*exp(-x.^2/2));
        for i = 0:order
            for j = 0:order
                basis = cat(3,basis,hg(i)'*hg(j));
            end
        end
    case 'cosine'
        for i = 0:order
            for j = 0:order
                basis = cat(3,basis,cos(pi*i*(0:window-1)/(window-1))'*cos(pi*j*(0:window-1)/(window-1)));
            end
        end
    otherwise
        error('Perhaps I will try some other basis functions at another time, but Hermite-Gaussian modes work pretty well for now')
end