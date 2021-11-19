function [x, output]=hiosd_lobpcg_sieh(F, H, x, options)
% semi-impicit HiOSD
global N c ep al epmckt2 dof kt
%% reading parameters
n=length(x);
if isfield(options,'innpfunc')
    innp = @(x,y)options.innpfunc(x,y);
    nor = @(x)sqrt(innp(x,x));
else
    disp('No inner product.');
    return
end
if isfield(options,'precd')
    P= options.precd;
else
    P=@(x)x;
end
options.k = size(options.V,2);
k = options.k;
V = options.V;
for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end
global X
dt=options.dt;
epsf=options.epsf;
maxiter=options.maxiter;

%% setting environment
alph=zeros ( k, 1) ;
U = zeros(n,k);
W = zeros(n,k);
Y = zeros(n,k);
%% bb iterations
Vp = zeros(n,k);
Up = zeros(n,k);
f = F (x);
gp = zeros(n,1);
Dx = f - 2 * V * innp(V,f); Dx = Dx * dt;
for iter=1:maxiter
    emp = true( 3 * k, 1 );
    tmpp = 2 * V * innp(V,f);
    g = f - tmpp;
    Dg = g - gp;
    gp = g;
    bta = abs ( innp(Dx,Dg) / innp(Dg,Dg) );
    bta = min ( bta, options.betat*dt);
    bta = max ( bta, options.betau*dt);
    
    xp = x;
    phi1 = [x; -sum(x)];
    phi2 = phi1 .^ 2 ;
    phi3 = phi2 .* phi1;
    ut = ifft2(reshape(phi1, N, N));  ut(1) = 0;
    
    nln = al * phi2 - phi3;
    nln = reshape(nln - mean(nln), N, N);
    
    ut = (ut + bta * ifft2(nln) - bta * ifft2(reshape([tmpp; -sum(tmpp)],N,N)) )...
        ./( 1 - bta * epmckt2 );    % semi-implicit
    
    x = real(fft2(ut)); x = x(1:dof)';
    Dx = x - xp;
    
    for i = 1 : k
        U(:,i) = H (x, V(:,i));
        alph(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alph(i) * V(:,i) ) ;
        W(:,i) = P(W(:,i)); 
        W(:,i) = W(:,i) - V * innp( V, W(:,i) ) ;
        W(:,i) = W(:,i) - W(:,1:i-1) * innp( W(:,1:i-1), W(:,i) ) ;
        nrmW = nor(W(:,i)) ;
        if nrmW > 1e-8
            W(:,i) = W(:,i) / nrmW ;
            Y(:,i) = H (x, W(:,i));
        else
            W(:,i) = 0;
            Y(:,i) = 0;
            emp(k+i) = false;
        end
    end
    for i = 1 : k
        zeta = innp( [V W Vp(:,1:i-1)], Vp(:,i) );
        Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * zeta;
        nrmVpi = nor(Vp(:,i));
        if nrmVpi > 1e-8
            Vp(:,i) = Vp(:,i) / nrmVpi ;
            Vp(:,i) = Vp(:,i) - [V W Vp(:,1:i-1)] * innp( [V W Vp(:,1:i-1)], Vp(:,i) );
            Vp(:,i) = Vp(:,i) / nor(Vp(:,i)) ;
            Up(:,i) = H (x, Vp(:,i));
        else
            Vp(:,i) = 0;
            Up(:,i) = 0;
            emp(2*k+i) = false;
        end
    end
    UU = [V,W,Vp]; UU = UU (:,emp);
    YY = [U,Y,Up]; YY = YY (:,emp);
    
    Pn = innp(UU,YY);  Pn = (Pn+Pn') / 2;
    [eta, alph] = eig( Pn, 'vector' );
    
    Vp=V;
    V = UU * real(eta(:,1:k));
    alph = real(alph(1:k));
    
    f = F (x);
    res = norm(ut + ifft2(nln) ./ epmckt2, 'fro') / N; 
    
    if res < epsf
        break;
    end
    if mod(iter,options.outputp)==0
        X=x;
        disp([num2str(iter) '   ' num2str(res) '   ' num2str(alph')]);
    end
    if mod(iter,options.outputd)==0
        eval(options.draws);
    end
end
output = struct('x', x, 'V', V, 'it', iter, 'alph', alph);

end
