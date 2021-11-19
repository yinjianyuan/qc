function [phi] = gradientflow(phi, maxiter, cname)
if nargin<2
    maxiter = 1e3;
end
global N c ep al epmckt2 dt dof kt
tol = 1e-12;  res = Inf;  iter = 0;

while res >= tol && iter < maxiter
    iter = iter + 1;
    phi1 = [phi; -sum(phi)];
    phi2 = phi1 .^ 2 ;
    phi3 = phi2 .* phi1;
    ut = ifft2(reshape(phi1, N, N));  ut(1) = 0;
    
    nln = al * phi2 - phi3;
    nln = reshape(nln - mean(nln), N, N);
    if mod(iter, 100) == 0
        res = norm(ut + ifft2(nln) ./ epmckt2, 'fro') / N;
        ene = (c/2) * norm(kt.*ut, 'fro')^2 + ...
            (-(ep/2).*sum(phi2) - (al/3).*sum(phi3) + 0.25.*norm(phi2, 'fro')^2) / (N^2);
        drawcam(phi); 
        fprintf('%d: %e, E = %.12e\n', iter, res, ene );
    end
    
    ut = (ut + dt*ifft2(nln) )./( 1 - dt*epmckt2 );    % semi-implicit
    phi1 = reshape(real(fft2(ut)), N^2, 1);
    phi = phi1(1:dof);
end
phi = phi1(1:dof);
if nargin == 3
    drawcam(phi, cname);
else
end

end