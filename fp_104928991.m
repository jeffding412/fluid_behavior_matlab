% Final Project
% Study fluid behavior in 2D using smoothed particle hydrodynamics (SPH) on
% a spatially-hashed domain
% Zhengfu Ding 104928991

% clear workspace
clear all;
clc;

N = 100; % number of particles
 
% create N particles 
particles(1:N) = struct('x',[],'y',[],'vel',struct('v_x',0,'v_y',0),'force',struct('f_x',0,'f_y',0),'rho',[],'neigh',[]);
 
% assign positions for each particle
for index = 1:N 
    particles(index).x = rand();
    particles(index).y = rand()/2;
end

% set the properties of the fluid being simulated
p0 = 1000; % Rest density
m = p0/N; % Mass of single particle
stiffness = 100; % Stiffness constant
mu = 0.1; % Viscosity coefficient
h = 0.25; % smoothing radius

% sets grid
xMax = 1; 
yMax = 1;
 
% bins
Nx = floor(xMax/h); % number of columns
Ny = floor(yMax/h); % number of rows
% bin dimensions
dx = xMax/Nx;
dy = yMax/Ny;

% construct Nx*Ny bins
numBins = Nx*Ny;
bins(1:numBins) = struct('particleIDs',[],'adjacentBins',[]);
 
% calculate the bin number for every particle
for k = 1:N
    binNum = (ceil(particles(k).x/dx)-1)*Ny + ceil((yMax-particles(k).y)/dy);
    bins(binNum).particleIDs = [bins(binNum).particleIDs, k];
end

% populate adjacentBins
for z = 1:numBins
    if z == 1 % top left corner
        bins(z).adjacentBins = [z+1,z+Ny,z+Ny+1];
    elseif z == Ny % bottom left corner 
        bins(z).adjacentBins = [z-1,z+Ny,z+Ny-1];
    elseif z == numBins-Ny+1 % top right corner
        bins(z).adjacentBins = [z+1,z-Ny,z-Ny+1];
    elseif z == numBins % bottom right corner
        bins(z).adjacentBins = [z-1,z-Ny,z-Ny-1];
    elseif z < Ny % left edge
        bins(z).adjacentBins = [z-1,z+1,z+Ny-1,z+Ny,z+Ny+1];
    elseif rem(z,Ny) == 1 % upper edge
        bins(z).adjacentBins = [z+1,z-Ny,z-Ny+1,z+Ny,z+Ny+1];
    elseif rem(z,Ny) == 0 % bottom edge
        bins(z).adjacentBins = [z-1,z-Ny,z-Ny-1,z+Ny,z+Ny-1];
    elseif z > numBins-Ny+1 % right edge
        bins(z).adjacentBins = [z+1,z-1,z-Ny,z-Ny+1,z-Ny-1];
    else % middle bins
        bins(z).adjacentBins = [z+1,z-1,z+Ny,z+Ny+1,z+Ny-1,z-Ny,z-Ny-1,z-Ny+1];
    end
end

% timestep
t = 1;
dt = 0.01;
damping = 0.75; % for wall collisions
 
% movie variables
vid = VideoWriter('hydrodynamics', 'MPEG-4');
vid.FrameRate = 30;
vid.Quality = 100;
% start video
open(vid);

while t < 5
    particles = identifyNeighbors(numBins, bins, h, particles, N);
    particles = calculateDensity(m, h, particles, N);
    particles = calculateForce(m, h, mu, particles, N);
    particles = updateKinematics(dt, particles, N, xMax, yMax, damping);
    
    % hold axes
    hold on;
    cla;
 
    % display data
    axis([0 xMax 0 yMax]);
    scatter([particles.x],[particles.y]);
    drawnow;
 
    % write to video
    writeVideo(vid, getframe(gcf));
    
    % timestep
    t = t + dt;
end

% end video
close(vid);

