clc; clear all; close all;
%% Script to render the perturbation trajectories of vortex system
% Read in vortex trajectory data
dir = '../DATA/';
fname = sprintf('%svortex_trajectories-GE-0500-005.x',dir);
fid = fopen(fname,'r','ieee-le');
nv = fread(fid,1,'int');
nt = fread(fid,1,'int');

yz  = zeros(2*nv,nt);
tau = zeros(1,nt);

for n = 1:nt
	yz(:,n) = fread(fid,2*nv,'double');
end

tau = fread(fid,nt,'double');

fprintf("Read in trajectory data\n");

%% Read in perturbation data
fname2 = sprintf('%sperturbations_3-0500-005.x',dir);
fid2 = fopen(fname2,'r','ieee-le');

nv = fread(fid2,1,'int');
nt = fread(fid2,1,'int');
nk = fread(fid2,1,'int');
tau = zeros(1,nt);
eta = zeros(nv,nt);
zeta = zeros(nv,nt);

yz_perturb = zeros(2*nv,nt);

tau = fread(fid2,nt,'double');

for n = 1:nt
	yz_perturb(:,n) = fread(fid2,2*nv,'double');
end

kb = fread(fid2,nk,'double');

eta(:,:)  = yz_perturb(1:nv,:);
zeta(:,:) = yz_perturb(nv+1:2*nv,:);

fprintf("Read in perturbation data\n");

%% Calculate theta for each vortex
theta = zeros(nv,nt);
for n = 1:nt
    for k = 1:nv
        theta(k,n) = atan(zeta(k,n)/eta(k,n));
    end
end

theta(2,:) = theta(2,:) - pi;
fprintf("Calculated theta displacement for each vortex\n");

%% Test section for plotting 3d sinusoidal line and rotating in global frame
np = 1000;
k = 0.73;
x = linspace(0, 2*pi/k,np);

y = zeros(nv,np,nt);
z = zeros(nv,np,nt);

Yp = zeros(1,np);
Zp = zeros(1,np);

amp = zeros(nv,nt);

for n = 1:nt
    for k = 1:nv
        amp(k,n) = 0.01*sqrt(eta(k,n)^2 + zeta(k,n)^2);
    end
end
%theta(:,:) = -theta(:,:);
for n = 1:nt
    for i = 1:nv
        Yp(:) = amp(i,n)*sin((k)*x);
        Zp(:) = 0.0;
        for j = 1:np
            [Yp(j), Zp(j)] = rotate_curve(Yp(j), Zp(j), theta(i,n));
            y(i,j,n) = yz(i,n) + Yp(j);
            z(i,j,n) = yz(i+nv,n) + Zp(j);
        end
    end
end

%% plot
figure(1)
for t_ind = 1:10:1500
clf;
plot3(x, y(1,:,t_ind), z(1,:,t_ind),'ro'), hold on
plot3(x, y(2,:,t_ind), z(2,:,t_ind),'bo'), hold on
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 2*pi/k])
ylim([-3 3])
zlim([0 5.5])
grid on
view(90,0)
%view(85,15)
m = 0.10;
pause(m); 
end

%% Test
t_ind = 900;

figure(2)
plot3(x, y(1,:,t_ind), z(1,:,t_ind),'ro'), hold on
plot3(x, y(2,:,t_ind), z(2,:,t_ind),'ko'), hold on
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 2*pi/k])
ylim([-3 3])
zlim([0 5.25])
grid on
view(90,0)
%view(85,15)



% 
% figure(2)
% for j = 1:10:nt
% t_ind = j;
% clf;
% plot3(x, y(t_ind,:,1), z(t_ind,:,1),'ro'), hold on
% plot3(x, y(t_ind,:,2), z(t_ind,:,2),'bo'), hold on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% xlim([0 2*pi/k])
% ylim([-2 2])
% zlim([0 5.25])
% grid on
% %view(90,0)
% view(80,15)
% m = 0.0001;
% pause(m);
% 
% end

%%
figure(3)
plot(tau,eta(1,:),'r-o'), hold on
plot(tau,zeta(1,:),'k-o'), hold on
plot(tau,eta(2,:),'r-'), hold on
plot(tau,zeta(2,:),'b-')
ylim([-30 30])
xlim([0 6.5])