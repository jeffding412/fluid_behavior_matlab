function [particles] = identifyNeighbors(numBins, bins, h, particles, N)
% populates each particle's neigh data member
% 5 inputs: 
% a scalar, numBins, representing the number of bins
% an array, bins, representing all bins
% a scalar, h, representing smoothing radius
% an array, particles, representing each individual particle
% and a scalar, N, representing number of particles
% returns array, particles, with neigh data
% Zhengfu Ding 104928991

% identify neighbors
for z = 1:numBins
    particlesInZ = bins(z).particleIDs; % get all particles in bin z
    adjacentBinsZ = bins(z).adjacentBins; % get vector of bin IDs adjacent to bin Z
    lengthAdjacentBins = length(adjacentBinsZ);
    for w = 1:lengthAdjacentBins
        particlesInW = bins(adjacentBinsZ(w)).particleIDs; % get all particles in bin W
        closeParticles = [particlesInW, particlesInZ]; % combine particles in bin W and bin Z
        for k = 1:length(particlesInZ) % all particles in bin Z
            for j = 1:length(closeParticles) % all particles in bin W
                % calculate distance
                dist = sqrt((particles(closeParticles(j)).x-particles(particlesInZ(k)).x)^2 + (particles(closeParticles(j)).y-particles(particlesInZ(k)).y)^2);
                if dist < h && closeParticles(j) ~= particlesInZ(k)
                    % particles k and j are neighbors
                    particleK = particlesInZ(k);
                    particleJ = closeParticles(j);
                    particles(particleK).neigh = [particles(particleK).neigh, particleJ];
                end
            end
        end
    end
end

% eliminate duplicate neighbors
for index = 1:N
    particles(index).neigh = unique(particles(index).neigh);
end
