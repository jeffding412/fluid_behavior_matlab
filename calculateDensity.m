function [particles] = calculateDensity(m, h, particles, N)
% calculates density of particles
% 4 inputs: 
% a scalar, m, representing the mass of a single particle
% a scalar, h, representing smoothing radius
% an array, particles, representing each individual particle
% and a scalar, N, representing number of particles
% returns array, particles, with density
% Zhengfu Ding 104928991

% calculate density
first = (4*m)/(pi*h*h);
secondCo = (4*m)/(pi*h^8);
for index = 1:N
    currentParticle = particles(index);
    neighborsLength = length(currentParticle.neigh);
    secondSum = 0;
    % iterate through neighbors
    for neighborIndex = 1:neighborsLength
        currentNeighbor = particles(currentParticle.neigh(neighborIndex));
        secondSum = secondSum + (h^2 - (sqrt((currentParticle.x-currentNeighbor.x)^2 + (currentParticle.y-currentNeighbor.y)^2))^2)^3;
    end
    particles(index).rho = first + secondCo*secondSum;
end
