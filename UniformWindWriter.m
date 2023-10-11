%% This file writes the uniform wind flow file for OpenFAST 
% This code takes in previously generated velocity data from TurbSim 2.0
% and calculates the Rotor-Effective Wind Speed (REWS). It then outputs
% REWS to a new wind file that can be used into OpenFAST. Be sure to change
% the variables WindType -> 2 and Filename_Uni -> "windfile.dat" inside
% InflowWind.dat 

% Written by Michael Lingad on 10/11/23

clc
clear
close all

%% Get the velocity and turbine data

[velocity, twrVelocity, y, z, zTwr, nz, ny, dz, dy, dt, zHub, z1,mffws] = readfile_BTS('90m_12mps_twr.bts');
r0 = 63;
HubCood = [0 90];
veclength = size(velocity,1);

%% Calculate Rotor-Effective Wind Speed (REWS)
Circle_angle = 0:pi/180:2*pi;
Circle_y = HubCood(1) + r0*cos(Circle_angle)' ; % Hub + R*cos(theta)
Circle_z = HubCood(2) + r0*sin(Circle_angle)' ; % Hub + R*sin(theta)
Circle_in = zeros(length(y),length(z));

% Determine the points within the rotor
for indy = 1:length(y)
    ypoint = y(indy);
    for indz = 1:length(z)
        zpoint = z(indz);
        Circle_in(indy,indz) = inpolygon(ypoint,zpoint,Circle_y,Circle_z);
    end
end

% Calculate REWS
REWS = zeros(1,veclength);
for ind = 1:veclength
    timeslice1 = squeeze(velocity(ind,1,:,:));
    rotortime = Circle_in.*timeslice1;
    rotortime = nonzeros(rotortime(:));
    REWS(ind) = sum(rotortime)/length(rotortime);
end

%%  Set up the vectors for the table
time = (-6.05:0.05:156)'; % Matching the time range found in ...\SFunc.IfW.Sum
windSpeed = REWS';
windDir = zeros(veclength,1);
vertSpeed = zeros(veclength,1);
horizShear = zeros(veclength,1);
pwrLawVertShr = zeros(veclength,1);
linVertShear = zeros(veclength,1);
gustSpeed = zeros(veclength,1);
upflowAngle = zeros(veclength,1);


% Define the output file name
outputFileName = 'REWSwind.dat';

% Open the file for writing
fid = fopen(outputFileName, 'w');

% Write the header
fprintf(fid, '! OpenFAST Deterministic Wind File\n');
fprintf(fid, '#\n');
fprintf(fid, '# Comment lines begin with "!" or "#" or "%%", then the data lines must contain the following columns:\n');
fprintf(fid, '#\n');
fprintf(fid, '# If there are only 8 columns, upflow is assumed to be 0.\n');
fprintf(fid, '#\n');
fprintf(fid, '# Parameters are interpolated linearly between time steps; using nearest neighbor before the first time\n');
fprintf(fid, '# listed in this file and after the last time listed in the file.\n');
fprintf(fid, '! Time     Wind    Wind    Vertical    Horiz.      Pwr.Law     Lin.Vert.   Gust     Upflow\n');
fprintf(fid, '!          Speed   Dir     Speed       Shear       Vert.Shr    Shear       Speed    Angle \n');
fprintf(fid, '! (sec)    (m/s)   (Deg)   (m/s)                                            (m/s)   (deg)\n');

% Combine the data into a matrix
data = [time, windSpeed, windDir, vertSpeed, horizShear, pwrLawVertShr, linVertShear, gustSpeed, upflowAngle];

% Write the data
for i = 1:size(data, 1)
    fprintf(fid, '%6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n', data(i,:));
end

% Close the file
fclose(fid);

disp(['Data has been written to ' outputFileName]);
