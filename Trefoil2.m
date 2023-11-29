% Define the parameter range
t = linspace(0, 100*pi, 10000);

% Define the height of the trefoil cone
h = 5;

% Parametric equations for a trefoil knot in x and y directions
x_knot = (sin(t) + 2 * sin(2*t));
y_knot = (cos(t) - 2 * cos(2*t));

% Define a scaling factor for the z-axis (you can modify this as needed)
scaling_factor = linspace(0, 1, length(t));

% Calculate the trefoil cone coordinates
x = x_knot .* scaling_factor;
y = y_knot .* scaling_factor;
z = h * t / max(t);
trefoil_data = [x', y', z'];

% Plot the trefoil cone
figure;
plot3(x, y, z, 'b', 'LineWidth', 2);
title('Trefoil Cone with Spiral Structure');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
grid on;
axis equal;

%% Based on Rodrigues rotation matrix
% R=I+sin(Δt⋅∥ω∥)⋅K+(1−cos(Δt⋅∥ω∥))⋅K2
% Assume 'trajectory' is an Nx3 matrix containing x, y, z coordinates at each point
N = size(trefoil_data, 1);

% Initialize cell arrays to store orientation matrices and angular velocities
orientation_matrices = cell(N, 1);
angular_velocities = cell(N, 1);

% Assume a small time step, you might need to adjust this based on your data
delta_t = 0.01;

% Iterate through the trajectory
for i = 1:N
    % Calculate the derivative of the position vector (velocity)
    if i == 1
        velocity = (trefoil_data(i+1, :) - trefoil_data(i, :)) / delta_t;
    elseif i == N
        velocity = (trefoil_data(i, :) - trefoil_data(i-1, :)) / delta_t;
    else
        velocity = (trefoil_data(i+1, :) - trefoil_data(i-1, :)) / (2 * delta_t);
    end

    % Calculate angular velocity
    angular_velocity = cross(trefoil_data(i, :), velocity) / norm(cross(trefoil_data(i, :), velocity));

    % Store angular velocity for further analysis if needed
    angular_velocities{i} = angular_velocity;

    % Convert angular velocity to rotation matrix
    K = [0, -angular_velocity(3), angular_velocity(2);
         angular_velocity(3), 0, -angular_velocity(1);
         -angular_velocity(2), angular_velocity(1), 0];

    rotation_matrix = eye(3) + sin(delta_t * norm(angular_velocity)) * K + (1 - cos(delta_t * norm(angular_velocity))) * K^2;

    % Store rotation matrix
    orientation_matrices{i} = rotation_matrix;
end

% 'orientation_matrices' now contains the orientation matrices at each point on the trajectory
%% Transformation matrix at each point
% Initialize cell array to store transformation matrices
transformation_matrices = cell(N, 1);

% Iterate through the trajectory
for i = 1:N
    % Extract orientation matrix for the current point
    orientation_matrix = orientation_matrices{i};

    % Create a 4x4 transformation matrix by concatenating [1,0] to trefoil_data row
    transformation_matrix = eye(4);
    transformation_matrix(1:3, 1:3) = orientation_matrix;
    transformation_matrix(1:3, 4) = trefoil_data(i, :)';

    % Store the transformation matrix
    transformation_matrices{i} = transformation_matrix;    

end
% 'transformation_matrices' now contains the 4x4 transformation matrices at each point on the trajectory
figure
for i = 1:500   
    coordplot(transformation_matrices{i},0.015)
    axis square
end
hold on