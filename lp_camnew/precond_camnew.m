function w = precond_camnew(w)
global N ikt2 dof
w=reshape([w;-sum(w)],N,N);%N
w=ifft2(w);
w=w.*ikt2;
w=real(fft2(w));
w=w(1:dof)';

end

