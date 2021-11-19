E=@ene_cammew;
F=@ngrad_camnew;
H=@hv_camnew;

options = hiosdoptions(1);
options.dt = 0.3*dt;     
options.betat = 8;       
options.betau = 0.2;
options.innpfunc = inpfp;
options.precd = @precond_camnew;
options.epsf = 1e-8; 
options.maxiter = 1000000;
options.outputp = 100;
options.outputd = 1000;   
options.draws='drawcam(x);';
options=rmfield(options,{'gammamax', 'gammamin', 'l', 'minl', 'tau'});
