function [particles] = calculateForce(m, h, mu, particles, N)
% calculates forces on particles
% 5 inputs: 
% a scalar, m, representing the mass of a single particle
% a scalar, h, representing smoothing radius
% a scalar, mu, representing the viscosity coefficient
% an array, particles, representing each individual particle
% and a scalar, N, representing number of particles
% returns array, particles, with forces
% Zhengfu Ding 104928991

% calculate force
for index = 1:N
    currentParticle = particles(index);
    neighborsLength = length(currentParticle.neigh);
    external = [0,-9.8]*currentParticle.rho + [1,0];
    sum = [0,0];
    % iterate through neighbors
    for neighborIndex = 1:neighborsLength
        currentNeighbor = particles(currentParticle.neigh(neighborIndex));
        % pressure
        first = (-30*m)/(pi*currentNeighbor.rho*h^4);
        second = (currentParticle.rho+currentNeighbor.rho)/2;
        qkj = (sqrt((currentParticle.x-currentNeighbor.x)^2 + (currentParticle.y-currentNeighbor.y)^2))/h;
        third = ((1-qkj)^2)/qkj;
        fourth = [currentParticle.x - currentNeighbor.x, currentParticle.y - currentNeighbor.y];
        sum = sum + first*second*third*fourth;
        % viscosity
        first = (-40*m*mu)/(pi*currentNeighbor.rho*h^4);
        second = 1 - qkj;
        third = [currentParticle.vel.v_x - currentNeighbor.vel.v_x, currentParticle.vel.v_y - currentNeighbor.vel.v_y];
        sum = sum + first*second*third;
    end
    finalVector = external + sum;
    particles(index).force.f_x = finalVector(1);
    particles(index).force.f_y = finalVector(2);
end
