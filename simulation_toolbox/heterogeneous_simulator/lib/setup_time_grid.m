function [kgrid, dt, Nt, pointsPerPeriod] = setup_time_grid(kgrid, medium, Nx, Ny, Nz, dx, dy, dz, frequency, CFL, varargin)
%SETUP_TIME_GRID Set up the time grid for k-Wave simulation.
%   [kgrid, dt, Nt, pointsPerPeriod] = setup_time_grid(kgrid, medium, Nx, Ny, Nz, dx, dy, dz, frequency, CFL)
%   [kgrid, dt, Nt, pointsPerPeriod] = setup_time_grid(..., 'display', true)
%
%   This function calculates and sets up the time grid parameters for a k-Wave
%   simulation based on the grid dimensions, medium properties, frequency, and CFL number.
%
%   Inputs:
%       kgrid      - kWaveGrid object
%       medium     - Structure with sound_speed field (3D array)
%       Nx, Ny, Nz - Grid dimensions
%       dx, dy, dz - Grid spacing in meters
%       frequency  - Source frequency in Hz
%       CFL        - Courant-Friedrichs-Lewy number
%
%   Optional Name-Value Pair Arguments:
%       'display'  - Display PPW and CFL information (default: true)
%
%   Outputs:
%       kgrid          - Updated kWaveGrid object with time set
%       dt             - Time step in seconds
%       Nt             - Number of time steps
%       pointsPerPeriod - Number of points per period
%
%   Example:
%       [kgrid, dt, Nt, pointsPerPeriod] = setup_time_grid(kgrid, medium, Nx, Ny, Nz, ...
%           dx, dy, dz, frequency, CFL);

% Parse optional arguments
p = inputParser;
p.addParameter('display', true, @islogical);
p.parse(varargin{:});
params = p.Results;

% Calculate reference and maximum sound speeds
c_ref = mean(medium.sound_speed(:));
c_max = max(medium.sound_speed(:)); % in m/s

% Calculate maximum distance in the grid
lmax = sqrt(((Nx-1)*dx)^2 + ((Ny-1)*dy)^2 + ((Nz-1)*dz)^2);

% Calculate end time
tend = lmax/c_ref;

% Calculate points per period based on CFL condition
pointsPerPeriod = round(c_max*1/frequency/dx/CFL);

% Calculate time step
dt = 1 / (pointsPerPeriod * frequency);

% Create the time array using an integer number of points per period
Nt = round(tend / dt);

% Set time on kgrid
kgrid.setTime(Nt, dt);

% Display information if requested
if params.display
    % Calculate the actual CFL and PPW
    disp(['PPW@ref= ' num2str(c_ref / (dx * frequency))]);
    disp(['CFL = ' num2str(c_max * dt / dx)]);
end

end

