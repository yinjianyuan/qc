function [ene] = ene_cammew(phi)
global N c ep al kt
phi1 = [phi; -sum(phi)];
phi2 = phi1 .^ 2 ;
ut = ifft2(reshape(phi1, N, N));  ut(1) = 0;

ene = (c/2) * norm(kt.*ut, 'fro')^2 + ...
    (-(ep/2).*sum(phi2) - (al/3).*sum(phi2 .* phi1) + 0.25.*norm(phi2, 'fro')^2) / (N^2);
end