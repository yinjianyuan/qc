function F = ngrad_camnew(phi)
global epmckt2 al N dof
phi = [phi; -sum(phi)];  
F = reshape( real(fft2(epmckt2.*ifft2(reshape(phi,N,N)))),N^2,1)...
    + (phi .^ 2) .* (al - phi);
F = F - mean(F);
F = F(1:dof);
end

