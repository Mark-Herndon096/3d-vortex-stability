clc; clear all; close all;


dir = '/Users/markherndon/3d-vortex-stability/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%somega.x',dir);
fid = fopen(fname,'r','ieee-le');

nv  = fread(fid,1,'int');
nk  = fread(fid,1,'int');
kb  = zeros(nk,1);
omg = zeros(nv,nk);

a = 0.098;

kb = fread(fid,nk,'double');

for k = 1:nk
	omg(:,k) = fread(fid,nv,'double');
end

ka = zeros(nk,1);
ka = kb*a;

%%
figure(1)
plot(ka,omg(1,:),'r-','LineWidth',1.5)
xlim([0 3])
ylim([0 1])
% 
% for k = 1:nk
%     bk0(k,1) = besselk(0,wvs(k));
%     bk1(k,1) = besselk(1,wvs(k));
%     bk2(k,1) = besselk(2,wvs(k));
% end
% figure(2)
% plot(wvs,bk0,'r-'), hold on
% plot(wvs,bk1,'b-'), hold on
% plot(wvs,bk2,'k-')
% ylim([-0.5 10])
