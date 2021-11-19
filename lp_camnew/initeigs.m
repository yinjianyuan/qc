E=@ene_cammew;
F=@ngrad_camnew;
H=@hv_camnew;

optionsini = hiosdoptions(1);
optionsini.dt = 0;       
optionsini.innpfunc=inpfp;
optionsini.outputX = 0; 
optionsini.outputp = 10;  
optionsini.outputd = 0;
optionsini.maxiter = 21;  
optionsini.precd = @precond_camnew;
optionsini.k=  10  ;   
rng(0);
optionsini.V=randn(dof,optionsini.k);

[V, alpha] = hiosd_lobpcg_iniP_eh(@hv_camnew,x,optionsini); 
optionsini.V=V; [V, alpha] = hiosd_lobpcg_iniP_eh(@hv_camnew,x,optionsini);
optionsini.maxiter = optionsini.maxiter *2-1;  
optionsini.V=V; [V, alpha] = hiosd_lobpcg_iniP_eh(@hv_camnew,x,optionsini);
optionsini.V=V; [V, alpha] = hiosd_lobpcg_iniP_eh(@hv_camnew,x,optionsini);
