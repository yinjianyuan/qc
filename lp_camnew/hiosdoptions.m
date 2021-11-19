function options = hiosdoptions(x)
options=struct('k',1,'dt',0.1,'l',1e-6,'minl',1e-6,...
    'epsf',1e-10,'maxiter',1e3,...
    'betat',max(x,1),'betau',min(1,1/x),'gammamax',max(x,1),'gammamin',min(1,1/x),...
    'outputX',1, 'outputp',1, 'outputd', 1, 'draws', '',...
    'tau',0.5,'precd',@(p)p);

if x==Inf
    options.betat=Inf;
    options.betau=0.01;
    options.gammammax=Inf;
    options.gammammin=1;
end
end
