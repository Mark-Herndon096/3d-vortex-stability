clc; clear all; close all;
global L_val k_val;

%% Script to render the perturbation trajectories of vortex system
% Read in vortex trajectory data
t_ind = 4750;
dir = '../DATA/';
fname = sprintf('%svortex_trajectories-GE-0200-005.x',dir);
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
fname2 = sprintf('%sperturbations_3-GE-0200-005.x',dir);
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
        theta(k,n) = atan2(zeta(k,n),eta(k,n));
    end
end

%theta(2,:) = theta(2,:)
fprintf("Calculated theta displacement for each vortex\n");

%% Test section for plotting 3d sinusoidal line and rotating in global frame
np = 1000;
k = 0.73;
L = 5*pi;
k_val = k;
L_val = 1.5*L/k_val;
x = linspace(0, L/k,np);

y = zeros(nv,np,nt);
z = zeros(nv,np,nt);

Yp = zeros(1,np);
Zp = zeros(1,np);

amp = zeros(nv,nt);

for n = 1:nt
    for j = 1:nv
        amp(j,n) = 0.05*sqrt(eta(j,n)^2 + zeta(j,n)^2);
    end
end

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
L2 = 1.5*L;
[xx, yy] = meshgrid(-L2/k:0.1:L2/k); % Generate x and y data
zz = zeros(size(xx, 1)); % Generate z data
figure(1)
for t_ind = 1:10:nt
clf;
surf(xx, yy, zz), hold on
shading interp
colormap bone
plot3(x, y(1,:,t_ind), z(1,:,t_ind),'bo','LineWidth',1.5), hold on
plot3(x, y(2,:,t_ind), z(2,:,t_ind),'ro','LineWidth',1.5), hold on
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 L/k])
ylim([-2 2])
zlim([0 5])
grid on
%view(90,0)
view(100,9)
m = 0.001;
pause(m); 
end

%% Test


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
[xx, yy] = meshgrid(-3*pi/k:0.1:3*pi/k); % Generate x and y data
zz = zeros(size(xx, 1)); % Generate z data
surf(xx, yy, zz)
shading interp
colormap bone




% 

%%
% figure(3)
% plot(tau,eta(1,:),'r-o'), hold on
% plot(tau,zeta(1,:),'k-o'), hold on
% plot(tau,eta(2,:),'r-'), hold on
% plot(tau,zeta(2,:),'b-')
% ylim([-30 30])
% xlim([0 6.5])
figure(4)
plot(tau,rad2deg(theta(1,:))), hold on
plot(tau,rad2deg(theta(2,:)),'ro')

%% 
figure(1)
for t_ind = 1:50:nt
tubeplot([x; y(1,:,t_ind); z(1,:,t_ind)],0.05,32);
tubeplot([x; y(2,:,t_ind); z(2,:,t_ind)],0.05,32);
m = 0.0001;
pause(m);
clf;
end