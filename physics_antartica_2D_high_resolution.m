%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linear heat equation PDE solver
% by Georges Tod, 2019
% georges.tod@outlook.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc, clear
addpath(genpath('functions'))

%% loading some data
load('data/antartica_highRES.mat')
figure(1);clf
image(geometry,'CDataMapping','scaled')
resolution=1;    % infrared image resolution
dx=(resolution)*1E-3;

% grid & timestamps
xmax=size(geometry,1);
ymax=size(geometry,2);

t_star = [0:0.1:10]';
w_star = zeros(xmax*ymax,length(t_star));

data.timeEXP=t_star;

% mapping geometry
n=1;

for x=1:xmax
    for y=1:ymax
        if geometry(x,y) == 1
            data.model(x,y)=1;
            solid_coords(n,1)=x;
            solid_coords(n,2)=y;
            n=n+1;
        else
            data.model(x,y)=0;
        end
    end
end

%{l
figure(2);clf
colormap winter
image(geometry,'CDataMapping','scaled')
caxis([0 10])
colorbar
hold on
plot(solid_coords(:,2),solid_coords(:,1),'.')
%}

% domain to compute
n=1;
for x=1:xmax
    for y=1:ymax
        if x>1 && x<xmax && y>1 && y<ymax
        %if (x-round(xmax/2))^2+(y-round(ymax/2))^2<32
            part.coords(n,1)=x;
            part.coords(n,2)=y;
            
            n=n+1;
        end
    end
end
clear x, clear y

%% Parameters & simulation
beta(1)=10;         % conductivity (in W/m/K)
beta(2)=NaN;          % length in z for BC (mm)
beta(3)=-40E7;       % power per volume unit (W/m^3)
beta(4)=1;            % energy time exposure  (s)

simulation_parameters.dx=dx;
simulation_parameters.xmax=xmax;
simulation_parameters.ymax=ymax;

tic
[time,S]=simul_2D(beta, simulation_parameters,data);
toc

return

%% plotting exp vs. sim

figure(3);clf
colormap winter

for i=2%:length(time)
    Msim=S(:,:,i)-273.15; % from kelvins to degrees
    image(Msim,'CDataMapping','scaled')
    colorbar
    caxis([-10 10])
    title(sprintf('simulation, t = %0.1f s',time(i)))
    
    Msim2 = Msim';
    wstar_sim(:,i-1)=Msim2(:);
    t_star_sim(i-1,1) = time(i-1);
    
    pause(0.1)
    
    
end


%% generating some data for machine learning

% training data - 2/3
t_star=t_star_sim(1:end*2/3);
w_star=wstar_sim(:,1:end*2/3)+273.15;

n=1;
for x=1:xmax
    for y=1:ymax
        
        X_star(n,1)=x;
        X_star(n,2)=y;
        
        n=n+1;
        
        
    end
end
%
save('antartica2DHR_training.mat','t_star','w_star','X_star')



% validation data - +1/3
t_star=t_star_sim;
w_star=wstar_sim+273.15;
%
save('antartica2DHR_validation.mat','t_star','w_star','X_star')


return

%% verification

figure(4);clf
colormap winter

for i=2%:length(time)
    Msim=reshape(w_star(:,i),[ymax,xmax])
    image(Msim,'CDataMapping','scaled')
    colorbar
    %caxis([0 10])
    
    
    
    
    
end


