clc; clear all; close all;


dir = '/Users/markherndon/3d-vortex-stability/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%somega.x',dir);
fname2 = sprintf('%somega_10.x',dir);
fid = fopen(fname2,'r','ieee-le');
fid2 = fopen(fname2,'r','ieee-le');

nk  = fread(fid,1,'int');
wvs = zeros(nk,1);
omg = zeros(nk,1);

nk2  = fread(fid2,1,'int');
wvs2 = zeros(nk2,1);
omg2 = zeros(nk2,1);

a = 0.31;

wvs = fread(fid,nk,'double');
omg = fread(fid,nk,'double');

wvs2 = fread(fid2,nk2,'double');
omg2 = fread(fid2,nk2,'double');
%%
figure(1)
plot(wvs(1:1:end),omg(1:1:end),'ko','LineWidth',1.5), hold on
plot(wvs2,omg2,'r-','LineWidth',1.5)
