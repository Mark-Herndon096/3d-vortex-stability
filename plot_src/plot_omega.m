clc; clear all; close all;


dir = '/Users/markherndon/3d-vortex-stability/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%somega.x',dir);
fid = fopen(fname,'r','ieee-le');


nk  = fread(fid,1,'int');
wvs = zeros(nk,1);
omg = zeros(nk,1);



a = 0.31;

wvs = fread(fid,nk,'double');
omg = fread(fid,nk,'double');

%%
figure(1)
plot(wvs(1:1:end),omg(1:1:end),'ko','LineWidth',1.5), hold on
plot(wvs,omg,'r-','LineWidth',1.5)
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
