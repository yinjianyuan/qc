clc
clear
addpath lp_camnew
close all
global Q q0 q1
global c ep al dt
global N L
%% Liquid to QC
Q = 12; q0 = 1; q1 = 2 * cos(pi/Q);  
c = 1; ep = -0.01; al = 1;         dt = 2;   
% choose the correct domain size:
L = 30; N = 256;
% L = 82;  N = 512;
% L = 112;  N = 1024; 

initialize_cam;
[x, cname] = guesses(7);
[x] = gradientflow(x, 10000, cname); % liquid
drawcam(x);drawnow

initeigs 
giveop;

x0=x+0.01*V(:,1);
options.V=V(:,[1]);
[x0, output]=hiosd_lobpcg_sieh(F, H, x0, options);% upward search a 1-saddle; 

V=output.V;
x01=x0+1e-5*V(:,1);
x02=x0-1e-5*V(:,1);
[x01] = gradientflow(x01, 10000);
[x02] = gradientflow(x02, 10000); % find the MEP from liquid to QC


%% QC to C6
Q = 12; q0 = 1; q1 = 2 * cos(pi/Q); 
c = 1; ep =  0.05; al = 1;         dt = 1;   
L = 82;  N = 512;
% L = 112;  N = 1024;

initialize_cam;
[x, cname] = guesses(1);
[x] = gradientflow(x, 10000, cname); % QC
drawcam(x);drawnow

ikt2 = 1 ./ ( kt.^2 + 0.04); 
initeigs 
giveop;
options.epsf = 1e-6; 

%upward search
x1=x+0.01*V(:,5);
options.V=V(:,[1:5]);
[x1, output]=hiosd_lobpcg_sieh(F, H, x1, options);

%downward search
V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1:4]);
x2=x1+0.001*V(:,5);
[x2, output]=hiosd_lobpcg_sieh(F, H, x2, options);

V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1:3]);
x3=x2+0.001*V(:,4);
[x3, output]=hiosd_lobpcg_sieh(F, H, x3, options);

V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1:2]);
x4=x3+0.001*V(:,3);
[x4, output]=hiosd_lobpcg_sieh(F, H, x4, options);

V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1]);
x5=x4+0.001*V(:,2);
[x5, output]=hiosd_lobpcg_sieh(F, H, x5, options);

% calculate mep
V=output.V;
x51=x5+1e-5*V(:,1);
x52=x5-1e-5*V(:,1);
[x51] = gradientflow(x51, 10000);
[x52] = gradientflow(x52, 10000);% find the MEP from QC to LQ


x=x52;
if ene_cammew(x52)>ene_cammew(x51)
    x=x51;
end 

initeigs 
giveop;
options.epsf = 1e-6; 

%upward search
x6=x+0.01*V(:,4);
options.V=V(:,[1:4]);
[x6, output]=hiosd_lobpcg_sieh(F, H, x6, options);

%downward search
V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1:3]);
x7=x6+0.001*V(:,4);
[x7, output]=hiosd_lobpcg_sieh(F, H, x7, options);

V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1:2]);
x8=x7+0.001*V(:,3);
[x8, output]=hiosd_lobpcg_sieh(F, H, x8, options);

V=output.V;
giveop;
options.epsf = 1e-7; 
options.V=V(:,[1]);
x9=x8+0.001*V(:,2);
[x9, output]=hiosd_lobpcg_sieh(F, H, x9, options);

% calculate mep
V=output.V;
x91=x9+1e-5*V(:,1);
x92=x9-1e-5*V(:,1);
[x91] = gradientflow(x91, 10000);
[x92] = gradientflow(x92, 10000);% find the MEP from LQ to C6





