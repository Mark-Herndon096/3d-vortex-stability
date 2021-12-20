
%% Script to render the perturbation trajectories of vortex system
clc; clear all; close all;

dir = '../DATA/';
fname = sprintf('%svortex_trajectories-GE-2000-005.x',dir);
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

% figure(1)
% plot(yz(1,:),yz(1+nv,:),'r-'), hold on
% plot(yz(2,:),yz(2+nv,:),'r-'), hold on
% ylim([0 6])
% xlim([-3 3])
% legend('Free vortex pair','Ground constraint','Center position','Center position');


%% Test section for plotting 3d sinusoidal line and rotating in global frame
np = 1000;
k = 0.73;
x = linspace(0, 2*pi/k,np);

y = zeros(nt,np,nv);
z = zeros(nt,np,nv);

Yp = zeros(1,np);
Zp = zeros(1,np);


%y = Yc + Yp;
%z(:) = Zc + Zp;

for n = 1:nt
    amp(n) = 0.01*exp(0.8*tau(n));
end


dtheta = 0.001;
theta(1) = 0 - 0.83;
theta(2) = pi + 0.83;

for n = 1:nv
for j = 1:nt
    Yp = amp(j)*sin(k*x);
    Zp(:) = 0.0;
    for i = 1:np
        [Yp(i), Zp(i)] = rotate_curve(Yp(i), Zp(i), theta(n));
        y(j,i,n) = yz(n,j) + Yp(i);
        z(j,i,n) = yz(nv+n,j) + Zp(i);
    end
end
end


figure(2)
for j = 1:10:nt
t_ind = j;
clf;
plot3(x, y(t_ind,:,1), z(t_ind,:,1),'ro'), hold on
plot3(x, y(t_ind,:,2), z(t_ind,:,2),'bo'), hold on
xlabel('x')
ylabel('y')
zlabel('z')
xlim([0 2*pi/k])
ylim([-2 2])
zlim([0 5.25])
grid on
%view(90,0)
view(80,15)
m = 0.0001;
pause(m);

end

