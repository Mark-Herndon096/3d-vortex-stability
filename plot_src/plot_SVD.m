%% Script to plot the perturbation trajectories of vortex system
clc; clear all; close all;

dir = '../DATA/';
fname = sprintf('%sperturbations_2-GE-0500-005.x',dir);
fid = fopen(fname,'r','ieee-le');

nv = fread(fid,1,'int');
nt = fread(fid,1,'int');
nk = fread(fid,1,'int');

tau = zeros(1,nt);
s = zeros(nk,nt);
kb = zeros(nk,1);

tau = fread(fid,nt,'double');

for n = 1:nt
	s(:,n) = fread(fid,nk,'double');
end

kb = fread(fid,nk,'double');


%% Plot trajectories
t_ind = 900;
% for k = 1:nk
%     s(k,t_ind) = log(s(k,t_ind))/tau(t_ind);
% end
figure(1)
plot(kb,s(:,t_ind),'k-')
%ylim([0 6])
xlim([0 6])
%legend('Free vortex pair','Ground constraint','Center position','Center position');
%set(L,'Interpreter','LaTeX')

[max_val, max_ind] = max(s(:,t_ind));
fprintf("s = %15.6f at kb = %5.4f\n",max_val, kb(max_ind))
