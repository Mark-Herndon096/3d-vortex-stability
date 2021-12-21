%% Script to plot the perturbation trajectories of vortex system
clc; clear all; close all;
t_ind = 1000;
dir = '../DATA/';
fname = sprintf('%sperturbations_3-GE-1000-005.x',dir);
fid = fopen(fname,'r','ieee-le');

nv = fread(fid,1,'int');
nt = fread(fid,1,'int');
nk = fread(fid,1,'int');
tau = zeros(1,nt);
eta = zeros(nv,nt);
zeta = zeros(nv,nt);

yz_perturb = zeros(2*nv,nt);

tau = fread(fid,nt,'double');

for n = 1:nt
	yz_perturb(:,n) = fread(fid,2*nv,'double');
end

kb = fread(fid,nk,'double');

eta(:,:)  = yz_perturb(1:nv,:);
zeta(:,:) = yz_perturb(nv+1:2*nv,:);

% plot
figure(1)
plot(tau,eta(1,:),'r-o'), hold on
plot(tau,zeta(1,:),'ko-'), hold on
plot(tau,eta(2,:),'r-'), hold on
plot(tau,zeta(2,:),'k-'), hold on
xlabel('Tau')
ylabel('Component amplitude')


amp = sqrt(yz_perturb(1,t_ind)^2 + yz_perturb(2,t_ind)^2 + yz_perturb(3,t_ind)^2 + yz_perturb(4,t_ind)^2);
fprintf('amp = %15.6f\n',amp);
