function hv = hv_camnew(phi,v)
global epmckt2 al N dof
v   = [  v; -sum(  v)];
phi = [phi; -sum(phi)];
hv = -reshape(real(fftn(epmckt2.*ifftn(reshape(v,N,N)))), N^2,1)...
    + (3*phi - 2*al) .* phi .* v;
hv = hv - mean(hv);
hv = hv(1:dof);
end

