%% Writing binary turbsim files using turbsim data (for testing ofc ofc)

clear
close all
clc

%% Loading stuff
tic
% load("readfile_BTS_Stuff.mat") % Tester bts file
load("LESValues_Yaw0_Rotor_1936.mat") % LES Rotor stuff
toc

%% LES Data Post-Processing
tmp = 7; % We're saving a non-periodic wind profile

% First, the LES data is nondimensionalized to 0.8, so divide everything by
% 0.8 and multiply by the rated velocity for the 5MW turbine which is 11.4 m/s
mffws = 20;

% Swap these two variables
Ny = ny; 
Nz = nz;
ny = Nz; % Number of spanwise points
nz = Ny; % Number of wall normal points
yhold = ytrunc;
zhold = ztrunc;
ytrunc = zhold; % Swap the y and z 
ztrunc = yhold;

% nt = size(U_rotor7,3); % Number of timesteps
d0  = 218; % Rotor diameter (5MW: 126, 15MW: 218)
dz = d0*(ztrunc(2) - ztrunc(1)); 
dy = d0*(ytrunc(2) - ytrunc(1));
dt = 2; % LES data has a timestep of 2 seconds, we might need to interpolate extra data...
zHub = 151; % Hub height (5MW: 90, 15MW: 150)
z1 = ztrunc(1)*d0 + zHub;

% Choose default Vslope Voffset values, we might need to fiddle with these later
Vslope = [1000 1000 1000];
Voffset = [0 0 0];

time_original = 0:dt:dt*(timesteps-1);
Deltat = 0.1; % New interpolate timestep
time_new = 0:Deltat:dt*(timesteps-1);
stepcount = length(time_new);
nt = stepcount;
U_interp = zeros(ny,nz,stepcount);
V_interp = zeros(ny,nz,stepcount);
W_interp = zeros(ny,nz,stepcount);
P_interp = zeros(ny,nz,stepcount);

% Let's choose x = 7 to feed to OpenFAST, this is 2D before the turbine! 
% I want to avoid the turbine effects to just get a "pure" flow
for iy = 1:ny
    for iz = 1:nz
        U_interp(iy,iz,:) = interp1(time_original,squeeze(U_rotor8(iy,iz,:)),time_new); % U_LES is streamwise, U_TurbSim is streamwise
        V_interp(iy,iz,:) = interp1(time_original,squeeze(W_rotor8(iy,iz,:)),time_new); % V_LES is wall-nornal, V_TurbSim is spanwise (need to swap!)
        W_interp(iy,iz,:) = interp1(time_original,squeeze(V_rotor8(iy,iz,:)),time_new); % W_LES is spanwise, W_TurbSim is wall-normal (need to swap!)
        P_interp(iy,iz,:) = interp1(time_original,squeeze(P_rotor8(iy,iz,:)),time_new); % P_LES data
    end
end
% save('P_interp.mat','P_interp')

% Undo the normalization and redimensionalize!
velocityU = mffws*U_interp/0.8;
velocityV = mffws*V_interp/0.8;
velocityW = mffws*W_interp/0.8;
velocity = permute(cat(4,velocityU,velocityV,velocityW),[3 4 1 2]); % Format this into [t,vel,y,z]

ntwr = length(find(ztrunc <= -0.5)); % Number of tower points BELOW the rotor
% Write the tower velocities
twrInd = find(ztrunc <= 5.01,1,'last'); % Turbine is located at [x = 9, y = 5, z = 0]

% Yoink the tower velocities
velocityUtwr = squeeze(velocityU(twrInd,1:ntwr,:));
velocityVtwr = squeeze(velocityV(twrInd,1:ntwr,:));
velocityWtwr = squeeze(velocityW(twrInd,1:ntwr,:));
twrVelocity = permute(cat(3,velocityUtwr,velocityVtwr,velocityWtwr),[2 3 1]); % Format this into [t,vel,z]

%% Setting up the writing to binary file process
tic
btsFileName = 'writeBTS_LESExpanded_at8_20mps_IEA15.bts';
fID_BTS = fopen(btsFileName,'W'); % Be sure to use 'W' for writing as this speeds up write times by A LOT

fwrite(fID_BTS,tmp,'int16'); % ID, should be 7 in general, 8 if the wind profile is periodic (INT(2))
fwrite(fID_BTS,nz,'int32'); % NumGrid_Z, Number of grid points in the wall-normal direction (INT(4))
fwrite(fID_BTS,ny,'int32'); % NumGrid_Y, Number of grid points in the spanwise direction (INT(4))
fwrite(fID_BTS,ntwr,'int32'); % ntower, Number of tower points below the grid (INT(4))
fwrite(fID_BTS,nt,'int32'); % n_t, Number of timesteps (INT(4))

fwrite(fID_BTS,dz,'float32'); % dz, Distance between spanwise points (z-direction) (REAL(4))
fwrite(fID_BTS,dy,'float32'); % dy, Distance between wall-normal points (y-direction) (REAL(4))
fwrite(fID_BTS,Deltat,'float32'); % TimeStep, Time in seconds between each slice (REAL(4))
fwrite(fID_BTS,mffws,'float32'); % Uhub, Mean wind speed at hub height in m/s (REAL(4)) 
fwrite(fID_BTS,zHub,'float32'); % HubHt, Height of the nacelle in m (REAL(4))
fwrite(fID_BTS,z1,'float32'); % Zbottom, Height in meters at the bottom of the grid (REAL(4))

% Write the slope and offset, this allows for the redimensionalization of velocity
% This allows for the conversion from 4 byte reals to 2 byte integers
fwrite(fID_BTS,Vslope(1),'float32'); % U-component slope for scaling (REAL(4))
fwrite(fID_BTS,Voffset(1),'float32'); % U-Component offset for scaling (REAL(4))
fwrite(fID_BTS,Vslope(2),'float32'); % V-component slope for scaling (REAL(4))
fwrite(fID_BTS,Voffset(2),'float32'); % V-Component offset for scaling (REAL(4))
fwrite(fID_BTS,Vslope(3),'float32'); % W-component slope for scaling (REAL(4))
fwrite(fID_BTS,Voffset(3),'float32'); % W-Component offset for scaling (REAL(4))

desc_string = 'This file was written in MATLAB using writefile_BTS.m and LES Data'; % Description string to send to TurbSim
char_count = length(desc_string); % Character count of the string
char_ascii = int32(double(desc_string)); % Translate this to ASCII

fwrite(fID_BTS,char_count,'int32') % Number of characters in ASCII (INT(4))

for ind = 1:char_count
    fwrite(fID_BTS,char_ascii(ind),'int8'); % Writing in the ASCII integer representation of the character string (INT(1))
end

disp('Writing velocities now')
for it = 1:nt % To number of timesteps
    % Write the rotor velocity
    for iz = 1:nz % To number of z points
        for iy = 1:ny % To number of y points
            for i = 1:3 % i = 1 -> U, i = 2 -> V, i = 3 -> W
                currVel = int16(velocity(it,i,iy,iz)*Vslope(i) + Voffset(i)); % Need to nondimensionalize to an int according to TurbSim Appendix D (D-1)
                fwrite(fID_BTS,currVel,'int16');
            end
        end
    end

    % Write the tower velocity
    for izTwr = 1:ntwr % To number of z tower points
        for i = 1:3 % i = 1 -> U, i = 2 -> V, i = 3 -> W
            currTwrVel = int16(twrVelocity(it,i,izTwr)*Vslope(i) + Voffset(i)); % Need to nondimensionalize to an int according to TurbSim Appendix D (D-2)
            fwrite(fID_BTS,currTwrVel,'int16');
        end
    end

    % Status line
    if mod(it,100) == 0
        disp(['Timesteps = ',num2str(it),' out of ',num2str(nt)])
    end
end

% The .bts file doesn't actually write the y and z vectors, it figures it
% out based on ny, dy, nz, and dz. This assumes they're equally spaced!
y = (0:ny-1)*dy - dy*(ny-1)/2;
z = (0:nz-1)*dz + z1;
zTwr = z1 - (0:ntwr-1)*dz;

fclose(fID_BTS);
toc

rrr

%% Check to see if the file can be read back to me
% tic
% [velocity, twrVelocity, yTS, zTS, zTwr, nz, ny, dz, dy, dt, zHub, z1,mffws] = readfile_BTS('90m_18mps_twr_30min_Coh_4D_Sim2.bts');
% toc

tic
[velocityTest, twrVelocityTest, yTSTest, zTSTest, zTwrTest, nzTest, nyTest, dzTest, dyTest, dtTest, zHubTest, z1Test,mffwsTest] = readfile_BTS(btsFileName);
toc

%% Plot the various velocities 
figure; % U Velocity (Streamwise Component)
for ind = 1:1:3001
    pcolor(y,z,squeeze(velocityU(:,:,ind))'); shading interp; colorbar vert; clim([6 15])
    title(['U Velocity, Time = ',num2str(time_new(ind))]); xlabel('y (Spanwise) [m]'); ylabel('z (Wall-Normal) [m]')
    pause(0.01)
end
% 
% figure; % V Velocity (Spanwise Component)
% for ind = 1:300
%     pcolor(y,z,squeeze(velocityV(:,:,ind))'); shading interp; colorbar vert;
%     title(['V Velocity, Time = ',num2str(time_new(ind))]); xlabel('y (Spanwise) [m]'); ylabel('z (Wall-Normal) [m]')
%     pause(0.1)
% end
% 
% figure; % W Velocity (Wall-Normal Component)
% for ind = 1:300
%     pcolor(y,z,squeeze(velocityW(:,:,ind))'); shading interp; colorbar vert;
%     title(['W Velocity, Time = ',num2str(time_new(ind))]); xlabel('y (Spanwise) [m]'); ylabel('z (Wall-Normal) [m]')
%     pause(0.1)
% end
% 
% figure; % Pressure Component
% for ind = 1:300
%     pcolor(y,z,squeeze(velocityW(:,:,ind))'); shading interp; colorbar vert;
%     title(['Pressure, Time = ',num2str(time_new(ind))]); xlabel('y (Spanwise) [m]'); ylabel('z (Wall-Normal) [m]')
%     pause(0.1)
% end
