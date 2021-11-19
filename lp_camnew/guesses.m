function [phi, cname] = guesses(choice)
initheta = 0;
sbid={
    [1:12];             %% 1 QC 
    [2:2:6];            %% 2 C6
    [7:2:11];           %% 3 CS6
    [2:2:6 1 7 12];     %% 4 LQ
    [8:2:12 11 2 3];    %% 5 LSQ
    [2 3 8];            %% 6 T6
    [];                 %% 7 H
    [1 ];               %% 8 Lam
    };
name = {'QC', 'C6', 'CS6', 'LQ', 'LSQ', 'T6', 'H', 'Lam'};
cname= [name{choice}];
disp(cname);

global Q q0 q1
global L N dof
kindex = [q0*[cos(initheta+      2*pi/12*(0:5)') sin(initheta+      2*pi/12*(0:5)')];
          q1*[cos(initheta+pi/12+2*pi/12*(0:5)') sin(initheta+pi/12+2*pi/12*(0:5)')]];
kindex = round([kindex; -kindex] * L); 
kindex(kindex<0) = kindex(kindex<0) + N;   kindex = kindex + 1;
ut = zeros(N,N);  
ut( kindex([sbid{choice} sbid{choice} + Q],:) * [1;N] -N ) = 0.058;

phi1 = reshape(real(fft2(ut)), N^2, 1);
phi = phi1(1:dof);

end