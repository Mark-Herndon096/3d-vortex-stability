%% Script to plot the perturbation trajectories of vortex system
clc; clear all; close all;

dir = '../DATA/';
fname = sprintf('%svortex_trajectories-2000-005.x',dir);
fid = fopen(fname,'r','ieee-le');
nv = fread(fid,1,'int');
nt = fread(fid,1,'int');

yz  = zeros(2*nv,nt);
tau = zeros(1,nt);

for n = 1:nt
	yz(:,n) = fread(fid,2*nv,'double');
end

tau = fread(fid,nt,'double');

%% Plot trajectories

figure(1)
plot(yz(1,:),yz(1+nv,:),'r-'), hold on
plot(yz(2,:),yz(2+nv,:),'r-'), hold on
plot(yz_2(1,:),yz_2(1+nv,:),'k-'), hold on
plot(yz_2(2,:),yz_2(2+nv,:),'k-'), hold on
ylim([0 6])
xlim([-3 3])
legend('Free vortex pair','Ground constraint','Center position','Center position');
%set(L,'Interpreter','LaTeX')
