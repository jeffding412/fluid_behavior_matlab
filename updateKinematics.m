function [particles] = updateKinematics(dt, particles, N, xMax, yMax, damping)
% updates positions and velocities of particles
% 6 inputs: 
% a scalar, dt, representing time step
% an array, particles, representing each individual particle
% a scalar, N, representing number of particles
% a scalar, xMax, representing the right boundary
% a scalar, yMax, representing the top boundary
% and a scalar, damping, for wall collisions
% returns array, particles, with updated positions and velocities using
% semi-explicit euler's method
% Zhengfu Ding 104928991

for index = 1:N
    % update velocity
   particles(index).vel.v_x = particles(index).vel.v_x + (dt*particles(index).force.f_x)/particles(index).rho;
   particles(index).vel.v_y = particles(index).vel.v_y + (dt*particles(index).force.f_y)/particles(index).rho;
   
   % update position
   particles(index).x = particles(index).x + dt*particles(index).vel.v_x;
   particles(index).y = particles(index).y + dt*particles(index).vel.v_y;
   
   % check for boundary behavior
   if particles(index).x > xMax % right boundary
       particles(index).x = 2*xMax - particles(index).x;
       particles(index).vel.v_x = -damping*particles(index).vel.v_x;
   end
   if particles(index).x < 0 % left boundary
       particles(index).x = 0 - particles(index).x;
       particles(index).vel.v_x = -damping*particles(index).vel.v_x;
   end
   if particles(index).y > yMax % top boundary
       particles(index).y = 2*yMax - particles(index).y;
       particles(index).vel.v_y = -damping*particles(index).vel.v_y;
   end
   if particles(index).y < 0 % bottom boundary
       particles(index).y = 0 - particles(index).y;
       particles(index).vel.v_y = -damping*particles(index).vel.v_y;
   end
end
