function drawcam(phi, NM)
global N L c ep al kt epmckt2 xs ys
phi1 = [phi; -sum(phi)];
ut = ifft2(reshape(phi1, N, N)); ut(1) = 0; 
utabs = abs(ut(:)); [~, index] = sort(utabs,'descend');
sample = min(128, sum(utabs(:)>1e-8)); index = index(1:sample)';

% set(gcf,'Position',[100 100 1500 800]);
set(gcf,'color','w');
subplot(1,2,1);
M1=L*pi; x0=linspace(-M1, M1, N+1); x0(end)=[];
pcolor(x0, x0, reshape(phi1, N, N)');shading interp;
set(gca,'Position',[0.01 0 0.48 1]);
xlim([-M1 M1]); ylim([-M1 M1]); axis equal; axis off;
ca=caxis; caxis([min(-0.0001,ca(1))  max(0.0001,ca(2))]);


phi2 = phi1 .^ 2 ;  phi3 = phi2 .* phi1;
ene = (c/2) * norm(kt.*ut, 'fro')^2 + ...
    (-(ep/2).*sum(phi2) - (al/3).*sum(phi3) + norm(phi2, 'fro')^2 /4) / (N^2);
nln = al * phi2 - phi3;		 nln = reshape(nln - mean(nln), N, N);
res = norm(ut + ifft2(nln)./epmckt2, 'fro') / N;
text(M1, M1*1.05,[ 'Energy = ' num2str(ene,'%.8g ')...
    '     normF = ' num2str(res ,'%.2e  ')],...
    'horiz','center','color',[0.5 0 0],'fontsize',16);
text(0, -M1*1.05,[ 'min = ' num2str(min(phi1),'%.3g ')...
    '     max = ' num2str(max(phi1) ,'%.3g  ')],...
    'horiz','center','color',[0 0 0],'fontsize',16);

subplot(1,2,2);
M2=4;  surf([-M2,M2],[-M2,M2],zeros(2),'FaceAlpha',1); shading flat; hold on; colormap parula; view(2);
utabs(index(abs(xs(index))>=M2))=NaN; utabs(index(abs(ys(index))>=M2))=NaN;
scatter3(xs(index), ys(index), utabs(index), 256, utabs(index), '.');
plot3(0,0,10,'.','Color',[1 1 1],'MarkerSize',16);
set(gca,'Position',[0.51 0 0.48 1]);
xlim([-M2 M2]);ylim([-M2 M2]);axis equal; axis off;
ca=caxis;caxis([0 max(0.0001,ca(2))]);
text(0, -M2*1.05,[ 'max = ' num2str(max(utabs(:)) ,'%.3g  ')],...
    'horiz','center','color',[0 0 0],'fontsize',16);
drawnow limitrate nocallbacks;
if nargin > 1
    print(gcf,['./lp_camnew/data' num2str(c) '_' num2str(ep) '_' num2str(al) '/' num2str(L) '_' num2str(N)  '/P' num2str(L) '_' NM],'-dpng','-r96','-opengl');
    save(     ['./lp_camnew/data' num2str(c) '_' num2str(ep) '_' num2str(al) '/' num2str(L) '_' num2str(N)  '/S' num2str(L) '_' NM], 'phi', '-ascii', '-double');
end
end
