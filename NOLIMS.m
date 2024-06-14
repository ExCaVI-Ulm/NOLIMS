%Non-Linear MRI Simulator; Version 06.2024

clear all; close all;close all hidden;

disp('Initialization')
addpath(genpath(pwd)); %add current path and all subfolders to search path
curr_path = pwd;    %needs to be where 'NOLIMS.m' is saved

%constants
settings.general.gamma = 2.675221874411*10^2;             % in rad/s/T; gyromagnetic ratio of H-nuclei.
settings.general.u0=1.26e-6;                              %u_0.
settings.general.T = 293.15;                              %Temperature in K: 20°C
settings.general.k = 1.380649e-23;                        %Boltzmann constant
settings.general.hbar = 1.05457182e-34;                   %hbar, reduced Planck's constant

%Field Setting
settings.general.B0 = 3;                                                                        %T mean Field--> in simArbFields3D updated but before used for Suscept.calculations, should be accurate for that
settings.general.FreqField = settings.general.gamma / (2*pi)*10^-6 * settings.general.B0;       %MHz; mean Freq of Field --> in ArbFieldsImport_MoreFlex Updated

%General settings for included effects
settings.general.noise = false;                                                                 % true: noise is added to simulated signal
settings.signal.SNR = 33;                                                                       %SNR in dB

%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE: not all features are available/tested in this example, it just demonstrates
%how time dependence might be integrated (in this 1D example)0 if desired in the basic parts of
%the software
%%%%%%%%%%%%%%%%%%%%%%%%

settings.general.T2 = false;                                                                    %true: simulate T2-decay: time-vector should not start at t=0;
settings.general.Suscept = false;                                                               %true: additional magnetic field due to susceptibility differences as additional phase during encoding: Assumption: B0 in z-dir (last dimension of susceptibility distribution)
settings.general.bool_IVD = false;                                                              %true: use sinc dephasing algorithm for signal generation: ATTENTION: then, matrixsize_signal = matrixsize_reco!; approximation: only gradients of field in z-direction are taken into account
settings.general.bool_IVD_reco = false;                                                         %true: incorporate dephasing (modeled with sincs) into reconstruction of data
settings.general.BlochSim = true;                                                              %true: Bloch Simulation for RX pulse with possibly underlying SEMs or Vector-character of B0
settings.general.CoilSens = false;                                                              %if false: uniform coil sens
settings.general.LowFieldw = false;                                                             %if true: w(r) in signal equation
settings.general.T1Effect = false;                                                              %if true: T1 effects between encoding steps with TR are included
settings.general.B0MapImport = false;                                                           %if true: import B0 Fieldmap
settings.general.B0MapInReco = false;                                                           %if true: B0Map in Reconstruction used
settings.general.B1MapImport = false;                                                           %if true: load B1Map for excitation, be careful with the format / left/right circularized etc.

settings.general.RAMSavingReco = false;                                                         %if true: ART reconstruction -> reduces RAM and enables GPU reco if not too large matrixsizes
settings.general.gpuComp = false;                                                               %if true: Reconstruction GPU for ART
settings.general.RAM_StepsPhaseEnc = 1;                                                       %if < encoding steps: instead of having all phase encoding steps in one matrix, divide them into multiple blocks, parameter indicates stepsize

settings.reco.FOV = 0.2;                        %FOV in m

settings.signal.matrixsize_signal = 64;         %matrix size used for simulating MRI signal

settings.reco.matrixsize_reco = 64;             %reconstructed matrix size: if matrixsize_signal > matrixsize_reco -> simulation of intravoxel dephasing -> factor between two values indicates how many subvoxels are simulated

%"trajectory"
settings.trajectory.RF_delay = 0;%0.08*10^-3;      %s; delay between RF and ACQ: Echo Time
settings.trajectory.TR = 600*10^-3;             %s; delay between subsequent RF excitations: Repetition Time

settings.trajectory.type = 'Rect1D';       %optional parameter to simulate different things

%%%%% IMPORTANT
%
%       THE DOMINANT COMPONENT OF THE IMPORTED MAGNETIC FIELD MUST BE IN
%       Z-DIRECTION / LAST DIMENSION OF THE TXT FILE / MATRIX
%
%%%%

settings.trajectory.dK = 1/settings.reco.FOV;
settings.trajectory.pixel_width = settings.reco.FOV / settings.reco.matrixsize_reco;
settings.trajectory.k_FOV = 1 / settings.trajectory.pixel_width;

settings.trajectory.Gread = 1*10^-3;            %T/m; Strength of Readout Gradient: fixed instead of readout length -> u.a. determined by fat water shift at scanner;

%Nyquist sampling
settings.trajectory.Tread = settings.trajectory.k_FOV *2*pi / ( settings.general.gamma * settings.trajectory.Gread );   %readout length

settings.trajectory.TOversamp = 1;  %Oversampling factor
settings.trajectory.BW = settings.trajectory.TOversamp*settings.general.gamma * settings.trajectory.Gread * settings.reco.FOV / (2*pi);% sampling Bandwidth


switch settings.trajectory.type
    case 'ImportFiles'
        settings.CST.Nx_CST = 32;                   %matrix size of field data exported with CST
        settings.CST.Ny_CST = settings.CST.Nx_CST;
        settings.CST.Nz_CST = settings.CST.Nx_CST;
        settings.CST.SixRows = 1; %if 1: data is saved with real and imag parts seperately  -> 6 export fields
        settings.CST.path_txtFiles = [curr_path, '/FieldData/'];     %path of folder where .txt files are saved
        
        %determine number of txt-files in folder to get number of rotations/"phase-encoding"-steps
        tmp_folderinfo = dir([settings.CST.path_txtFiles '/*.txt']);        
        settings.trajectory.N_PhaseEnc = length(tmp_folderinfo);
        
        settings.trajectory.BW =settings.trajectory.BW; %either take values calculated from Nyquist Sampling / standard encoding or set new values
        settings.trajectory.Tread = 0.5*settings.trajectory.Tread;  %above: "cartesian" calculation but it's rather UTE-like -> half readout length
        settings.trajectory.Nsamples = ceil(settings.trajectory.Tread * settings.trajectory.BW); %number of samples
        
        clearvars tmp_folderinfo
    case 'Rect1D'
        settings.trajectory.N_PhaseEnc = 1;
        settings.trajectory.BW =settings.trajectory.BW; %either take values calculated from Nyquist Sampling / standard encoding or set new values
        settings.trajectory.Tread = settings.trajectory.Tread;  %above: "cartesian" calculation but it's rather UTE-like -> half readout length
        settings.trajectory.Nsamples = (settings.trajectory.Tread * settings.trajectory.BW); %number of samples
    otherwise
        warning('Unrecognized sampling pattern')
end

if settings.general.B0MapImport
    pathB0 = [curr_path, '/FieldData/B0Map.mat'];     %path of B0Map
    ImportB0 = load(pathB0);

    %interpolate to requested size
    Nx_B0=size(ImportB0.B0Map,1); Ny_B0=size(ImportB0.B0Map,2); Nz_B0=size(ImportB0.B0Map,3);
    x_B0 = ([0:Nx_B0-1]-Nx_B0/2+0.5);y_B0 = ([0:Ny_B0-1]-Ny_B0/2+0.5);z_B0 = ([0:Nz_B0-1]-Nz_B0/2+0.5);
    X_B0=linspace(x_B0(1), x_B0(Nx_B0), settings.signal.matrixsize_signal);
    Y_B0=linspace(y_B0(1), y_B0(Ny_B0), settings.signal.matrixsize_signal);
    Z_B0=linspace(z_B0(1), z_B0(Nz_B0), settings.signal.matrixsize_signal);

    X_B0_Reco=linspace(x_B0(1), x_B0(Nx_B0), settings.reco.matrixsize_reco);
    Y_B0_Reco=linspace(y_B0(1), y_B0(Ny_B0), settings.reco.matrixsize_reco);
    Z_B0_Reco=linspace(z_B0(1), z_B0(Nz_B0), settings.reco.matrixsize_reco);
    
    B0MapInterp =interp3(y_B0, x_B0, z_B0, ImportB0.B0Map, Y_B0', X_B0, Z_B0, 'linear');
    B0MapInterp_Reco =interp3(y_B0, x_B0, z_B0, ImportB0.B0Map, Y_B0_Reco', X_B0_Reco, Z_B0_Reco, 'linear');

    settings.signal.B0Map = B0MapInterp;
    settings.reco.B0Map = B0MapInterp_Reco;
end

if mod(settings.trajectory.N_PhaseEnc, settings.general.RAM_StepsPhaseEnc) ~= 0 || settings.general.RAM_StepsPhaseEnc >  settings.trajectory.N_PhaseEnc
    error('StepSize not valid!!')
end

if settings.general.bool_IVD
    if settings.signal.matrixsize_signal ~= settings.reco.matrixsize_reco
        error('Invalid settings: IVD Dephasing but  different sizes for signal and reco')
    end
end

%Bloch Simulation of TX pulse
settings.TX.pulse_type = 'block';%Alternatives: block, sinc, gaussian10

settings.TX.PulseLength = 30.4*10^-6;%s
switch settings.TX.pulse_type   %Pulse integrals etc. from De Graaf, Robin A. In vivo NMR spectroscopy: principles and techniques. John Wiley & Sons, 2019.
    case 'block'
        settings.TX.Rp = 1.37;
        settings.TX.pulse_integral = 1;
    case 'sinc'
        settings.TX.Rp = 5.95;
        settings.TX.pulse_integral = 0.1775;
    case 'gaussian10'
        settings.TX.Rp = 2.0;
        settings.TX.pulse_integral = 0.565;
    otherwise
        warning('Not implemented');
end

settings.TX.PulseBW = settings.TX.Rp/settings.TX.PulseLength; %Hz:  bandwidth of exc. pulse
settings.TX.FlipAngle = 90*pi/180;   %flip angle in rad
settings.TX.PulseAmpl = settings.TX.FlipAngle/(settings.general.gamma*settings.TX.PulseLength*settings.TX.pulse_integral);  %amplitude of RF pulse
settings.TX.ExcSensitivity=1;
settings.TX.AngleRF_Theta = 0; %rad; phase of excitation pulse

settings.TX.NumberSamplesBloch = 1000;%numbers of time steps for Bloch Simulation

% RX Coils
settings.CoilSens.Ncoil_segments =550;                              %number of segments the coil is divided into 
settings.CoilSens.Curr=1;                                           %current in the coil   
settings.CoilSens.NReceiveCoils=8;                                  %number of der receiver coils
settings.CoilSens.RadiusReceiveCoils = 0.1;                         %m, Radius of circular loops
settings.CoilSens.RadiusCoilArray = (settings.reco.FOV - 0.005)/2;  %m, Radius of Coil Array: Coils 0.5cm away from edge of FOV
settings.CoilSens.DistMaskWire = 0.001;                             %treshold of distance between wire and position where not to evaluate B1 / CoilSens

%Sample
settings.sample.type = 'Rect1D'; %alternatives: SheppLogan, Sphere, Cylinder, Import
switch settings.sample.type
    case 'SheppLogan'
        
    case 'Sphere'
        settings.sample.sph_shiftx = 2.5*10^-3;         %m
        settings.sample.sph_shifty = -8.52*10^-3;       %m
        settings.sample.sph_shiftz = -8.52*10^-3;       %m
        settings.sample.sph_radius = 0.0725;            %m
        settings.sample.sph_T1 = 0.218;                 %s
        settings.sample.sph_T2 = 0.1735;                %s
        settings.sample.sph_SusceptExt = 0.36*10^-6;    %susceptibility outside the phantom
        settings.sample.sph_SusceptInt = -9.05*10^-6;   %susceptibility inside the phantom
    case 'Cylinder'

        settings.sample.cyl_shiftx = -0.05;             %m
        settings.sample.cyl_shifty = 0.01;              %m
        settings.sample.cyl_shiftz = 0.01;              %m
        settings.sample.cyl_radius = 0.04;              %m
        settings.sample.cyl_height = 0.16;              %m
        settings.sample.cyl_T1 = 2.75;                  %s
        settings.sample.cyl_T2 = 2.05;                  %s
        settings.sample.cyl_SusceptExt = 0.36*10^-6;    %susceptibility outside the phantom
        settings.sample.cyl_SusceptInt = -9.05*10^-6;   %susceptibility inside the phantom
    case 'Import'
        settings.sample.path = [curr_path, '/Sample/']; %path where to find the files of a sample
    case 'Rect1D'

    otherwise
        warning('Not implemented');
end

if ~settings.general.CoilSens
    settings.CoilSens.NReceiveCoils = 1; %minimize number of unnecessary for loops and matrixsizes if uniform coil sensitivity is desired
end

if settings.general.B1MapImport
    pathB1 = [curr_path, '/FieldData/B1Map.mat'];     %path of B1Map size: settings.trajectory.N_PhaseEnc,3,time, matrix spatial, matrix spatial, matrix spatial
    ImportB1 = load(pathB1);

    %interpolate to requested size
    Nx_B1=size(ImportB1.B1Map,1); Ny_B1=size(ImportB1.B1Map,2); Nz_B1=size(ImportB1.B1Map,3);
    x_B1 = ([0:Nx_B1-1]-Nx_B1/2+0.5);y_B1 = ([0:Ny_B1-1]-Ny_B1/2+0.5);z_B1 = ([0:Nz_B1-1]-Nz_B1/2+0.5);
    X_B1=linspace(x_B1(1), x_B1(Nx_B1), settings.signal.matrixsize_signal);
    Y_B1=linspace(y_B1(1), y_B1(Ny_B1), settings.signal.matrixsize_signal);
    Z_B1=linspace(z_B1(1), z_B1(Nz_B1), settings.signal.matrixsize_signal);
    ImB1 = ImportB1.B1Map;
    B1MapInterp = zeros(size(ImB1,1), size(ImB1,2), settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);
    
    for ll = 1:settings.trajectory.N_PhaseEnc   %different RF for each readout?
        for j =1:3  %x,y,z
            for u = 1:settings.TX.NumberSamplesBloch
                B1MapInterp(ll, j,u,:,:,:) =interp3(y_B1, x_B1, z_B1, squeeze(ImB1(ll,j,u,:,:,:)), Y_B1', X_B1, Z_B1, 'linear');
            end
        end
    end
    settings.signal.B1Map = B1MapInterp;

else    
    %if not imported: "build a simulation B1Map"
    B1 = zeros(3, settings.TX.NumberSamplesBloch, settings.signal.matrixsize_signal);  %x,y,z in first dimension; 2nd dimension time, last dim: spatial

    %Rest (time dependence / oscillation / frequency) in Bloch Sim function
    settings.signal.B1MapNoImport = B1;
end
%% Coordinates and Sample

FOV = settings.reco.FOV;
matrixsize_signal = settings.signal.matrixsize_signal;
matrixsize_reco = settings.reco.matrixsize_reco;

%Coordinate System
%Field is assumed to be in z-direction
%%%%%%%%%%%%%
% î
% x       ________    
%       _|______  |
%     _|_______ | |
%    |        | | |
%    |        | |_|    z
%    |        |_|      ^
%    |________|       /
%
% -> y
% first dim: x: 1 =pos ->end neg
%second dim: y: 1= neg -> end pos
%third dim: z: 1=neg -> end pos

%Coordinate systems if matrix is shown with as()
%if x constant:       | -> z
%                   y v 
%
%if y constant:      x   
%                     î -> z
%
%if z constant:   x
%                   î -> y

%%%%%%%%%%%%%

x = fliplr(linspace(-FOV/2, FOV/2- settings.trajectory.pixel_width * matrixsize_reco/matrixsize_signal, matrixsize_signal));
y = fliplr(x);
z = fliplr(x);

[coord_x, coord_y, coord_z] = meshgrid(x,y,z);
coord_x = permute(coord_x, [2 1 3]);
coord_y = permute(coord_y, [2 1 3]);
coord_z = permute(coord_z, [2 1 3]);

x_reco = fliplr(linspace(-FOV/2, FOV/2- settings.trajectory.pixel_width, matrixsize_reco));
y_reco = fliplr(x_reco);
z_reco = fliplr(x_reco);

[coord_x_reco, coord_y_reco, coord_z_reco] = meshgrid(x_reco,y_reco,z_reco);
coord_x_reco = permute(coord_x_reco, [2 1 3]);
coord_y_reco = permute(coord_y_reco, [2 1 3]);
coord_z_reco = permute(coord_z_reco, [2 1 3]);

%% Sample
disp('Sample');tic
air_suscept = 0.36*10^-6;   %susceptibility of air

switch settings.sample.type
    case 'SheppLogan'

        shepplogan = (phantom('Modified Shepp-Logan', matrixsize_signal));
        shepplogan = repmat(shepplogan, [1,1,matrixsize_signal]);   %make it 3D

        sampleS.M0 = shepplogan;

        sample_straight = reshape(shepplogan, 1, []);
        sample_straight(2,:) =  reshape(coord_x, 1, []); %x
        sample_straight(3,:) =  reshape(coord_y, 1, []); %y
        sample_straight(4,:) =  reshape(coord_z, 1, []); %z

        %Susceptibility and T2
        
        Suscept = [  (-8*10^-6 )  .69   .92    0     0     0             %lipids
                (-9.05*10^-6  -(-8*10^-6 ))  .6624 .8740   0  -.0184   0              %water
                (0- (-9.05*10^-6 )) .1600 .4100 -.22    0     18                                      %air
                (0- (-9.05*10^-6 ))  .1100 .3100  .22    0    -18
                (-9.05*10^-6  -(-9.05*10^-6 ))  .2100 .2500   0    .35    0              %water
                (-9.05*10^-6  -(-9.05*10^-6 ))  .0460 .0460   0    .1     0
                (-9.05*10^-6  -(-9.05*10^-6 ))  .0460 .0460   0   -.1     0
                (-7.9*10^-6 - (-9.05*10^-6 ))  .0460 .0230 -.08  -.605   0               %deoxy blood
                (-7.9*10^-6 - (-9.05*10^-6 ))  .0230 .0230   0   -.606   0
                (-7.9*10^-6 - (-9.05*10^-6 ))  .0230 .0460  .06  -.605   0   ];


        T2d = [  70 .69   .92    0     0     0   
                100-70  .6624 .8740   0  -.0184   0
                500-100 .1600 .4100 -.22    0     18
                500-100  .1100 .3100  .22    0    -18
                120-100  .2100 .2500   0    .35    0
                60-100  .0460 .0460   0    .1     0
                120-100  .0460 .0460   0   -.1     0
                80-100  .0460 .0230 -.08  -.605   0 
                80-100  .0230 .0230   0   -.606   0
                80-100  .0230 .0460  .06  -.605   0   ];

        %T1 only for Bloch-Simulation of TX pulse
        T1d = [  450  .69   .92    0     0     0   
                950-450  .6624 .8740   0  -.0184   0
                2000-950  .1100 .3100  .22    0    -18
                2000-950 .1600 .4100 -.22    0     18
                600-950  .2100 .2500   0    .35    0
                1200-950  .0460 .0460   0    .1     0
                1200-950  .0460 .0460   0   -.1     0
                1200-950  .0460 .0230 -.08  -.605   0 
                1200-950  .0230 .0230   0   -.606   0
                1200-950  .0230 .0460  .06  -.605   0   ];
        %T2
        P2 = phantom(T2d,matrixsize_signal);
        P2(P2==0) = Inf;    %T2 = 0 is replaced by Inf
        sampleT2 = P2*10^-3;
        sampleT2 = repmat(sampleT2, [1,1,matrixsize_signal]);   %make it 3D
        sampleS.T2 = sampleT2;
        sample_straight(5, :) = reshape(sampleT2, 1, []); %T2

        %T1
        P1 = phantom(T1d,matrixsize_signal);
        P1(P1==0) = Inf;    %T1 = 0 is replaced by Inf
        sampleT1 = P1*10^-3;
        sampleT1 = repmat(sampleT1, [1,1,matrixsize_signal]);   %make it 3D
        sampleS.T1 = sampleT1;
        sample_straight(6, :) = reshape(sampleT1, 1, []); %T1

        %Susceptibility: makes only limited sense -> rather use "localised" off-resonance
        PH_sus = phantom(Suscept, matrixsize_signal);
        PH_sus(PH_sus == 0) = air_suscept;
        PH_sus = PH_sus - air_suscept;  %sign of susceptibility difference "yields switched x, y coordinates", 1ppm as theory said, -1ppm rotated by 90deg
        PH_sus_3D = repmat(PH_sus, [1 1 matrixsize_signal]);

        dB = settings.general.B0*calculateFieldShift(PH_sus_3D, [FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal]*10^3); %https://de.mathworks.com/matlabcentral/fileexchange/56680-mri-simulation-using-forecast-fourier-based-off-resonance-artifact-simulation-in-the-steady-state
        DB_straight = reshape(squeeze(dB),1,[]);
        sampleS.DB = squeeze(dB);

        
        %Localised Off-Resonance instead of weird Off-Res for SheppLogan being only "2D"

        f0map = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        Radius_f0 = 0.042;%m
        shift_x0 =0;
        shift_y0 =-0.05;
        shift_z0 =0;
        I=(coord_x-shift_x0).^2+(coord_y-shift_y0).^2+(coord_z-shift_z0).^2< Radius_f0^2;
        f0map(I)=2800;%Hz
        as(f0map)
        f0map_blur = imgaussfilt(f0map,2);
        as(f0map_blur)

        dB = 2*pi*f0map / settings.general.gamma;
        DB_straight = reshape(squeeze(dB),1,[]);
        sampleS.DB = dB;
    case 'Sphere'
        SpherePhantom = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        
        x_s = x;
        y_s = y;
        z_s = z;
        [coord_s_x, coord_s_y, coord_s_z] = meshgrid(x_s,y_s, z_s);
        
        I_sphPh=(coord_s_x-settings.sample.sph_shiftx).^2+(coord_s_y-settings.sample.sph_shifty).^2 +(coord_s_z-settings.sample.sph_shiftz).^2< settings.sample.sph_radius^2; %
        
        SpherePhantom(I_sphPh)=1;%M0        
        sampleS.M0 = squeeze(SpherePhantom);
        
        sampleS.T1 = settings.sample.sph_T1*squeeze(SpherePhantom);
        sampleS.T1(sampleS.T1==0) = Inf;    %T1 = 0 is replaced by Inf
        sampleS.T2 = settings.sample.sph_T2*squeeze(SpherePhantom);
        sampleS.T2(sampleS.T2==0) = Inf;    %T2 = 0 is replaced by Inf

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(coord_s_x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(coord_s_y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(coord_s_z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
        Susc_sph = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        Susc_sph(I_sphPh) = abs(settings.sample.sph_SusceptInt - settings.sample.sph_SusceptExt);
        
        dB = settings.general.B0*calculateFieldShift(Susc_sph, [FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal]*10^3);
        DB_straight = reshape(squeeze(dB),1,[]);
        
        sampleS.DB = squeeze(dB);
        clearvars coord_s_x coord_s_y coord_s_z Susc_sph SpherePhantom x_s y_s z_s I_sphPh
    case 'Cylinder'
        CylinderPhantom = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        bool_cyl = (and((coord_z-settings.sample.cyl_shiftz).^2 + (coord_y-settings.sample.cyl_shifty).^2 <= settings.sample.cyl_radius^2, and(coord_x>settings.sample.cyl_shiftx,coord_x<(settings.sample.cyl_height+settings.sample.cyl_shiftx))));  %&& Z>0 && Z<height_cylinder
        
        CylinderPhantom(bool_cyl)=1;%M0        
        sampleS.M0 = squeeze(CylinderPhantom);
        
        sampleS.T1 = settings.sample.cyl_T1*squeeze(CylinderPhantom);
        sampleS.T1(sampleS.T1==0) = Inf;    %T1 = 0 is replaced by Inf --> leads to wrong values for scaling of T1 but there is no Magnetization either
        sampleS.T2 = settings.sample.cyl_T2*squeeze(CylinderPhantom);
        sampleS.T2(sampleS.T2==0) = Inf;    %T2 = 0 is replaced by Inf

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(coord_x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(coord_y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(coord_z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
        Susc_cyl = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        Susc_cyl(bool_cyl) = abs(settings.sample.cyl_SusceptInt - settings.sample.cyl_SusceptExt);
        
        dB = settings.general.B0*calculateFieldShift(Susc_cyl, [FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal,FOV/settings.signal.matrixsize_signal]*10^3);
        DB_straight = reshape(squeeze(dB),1,[]);
        sampleS.DB = squeeze(dB);

    case 'Import'
        %examplary load JEMRIS Brain Phantom
        load([settings.sample.path,'MNIbrain.mat'], 'BRAIN');
        load([settings.sample.path,'MNIdeltaB.mat'], 'DB');

        %Code from JEMRIS

        %tissue parameters from JEMRIS, mrm.29009 (mean) @50mT, Bottomley NMR Relax in tissue
        %assuming: T2 is independent of field strength
        %        T1  T2 T2*[ms]  M0 CS[rad/sec]      Label
        tissue=[3702 329  158   1.00   0         ;  % 1 = CSF
                 326.7  83   69   0.86   0         ;  % 2 = GM
                 274.7  70   61   0.77   0         ;  % 3 = WM
                 146  70   58   1.00 7.9*2*pi    ;  % 4 = Fat (CS @ 53mT Tesla) ->suscept map
                 130  47   30   1.00   0         ;  % 5 = Muscle / Skin
                 450 329   58   1.00   0         ;  % 6 = Skin
                   0   0    0   0.00   0         ;  % 7 = Skull
                 826  83   69   0.86   0         ;  % 8 = Glial Matter
                 122  70   61   0.77   0         ;];% 9 = Meat

        %parameter maps
        PARAMS={'M0','T1','T2','T2S','DB'};
        fact=[1 1 1 1 1]; %if 1: activated; M0, T1, T2, T2s, CS
        INDEX =[4 1 2 3 5];
        for i=1:9
            for j=1:5
                if i==1,eval(['BrainSample.',PARAMS{j},'=zeros(size(BRAIN));']);end
                I   = find(BRAIN==i);
                ind = INDEX(j);
                eval(['BrainSample.',PARAMS{j},'(I)=fact(j)*tissue(i,ind);']);
            end
        end

        %add susceptibility issues
        %interpolate DB to twice the initial size
        Nx_suscep=size(DB,1); x_suscep=([0:Nx_suscep-1]-Nx_suscep/2+0.5);
        Ny_suscep=size(DB,2); y_suscep=([0:Ny_suscep-1]-Ny_suscep/2+0.5);
        Nz_suscep=size(DB,3); z_suscep=([0:Nz_suscep-1]-Nz_suscep/2+0.5);
        X_suscep=[x_suscep(1):(x_suscep(Nx_suscep)-x_suscep(1))/(2*Nx_suscep-1):x_suscep(Nx_suscep)];
        Y_suscep=[y_suscep(1):(y_suscep(Ny_suscep)-y_suscep(1))/(2*Ny_suscep-1):y_suscep(Ny_suscep)];
        Z_suscep=[z_suscep(1):(z_suscep(Nz_suscep)-z_suscep(1))/(2*Nz_suscep-1):z_suscep(Nz_suscep)];

        DB = interp3(y_suscep,x_suscep,z_suscep,DB,Y_suscep',X_suscep,Z_suscep,'spline');

        BrainSample.DB = BrainSample.DB + 2*pi*1e6*DB*settings.general.FreqField;
        for j =1:5
            eval(['BrainSample.',PARAMS{j},'=flip(permute(shiftdim(BrainSample.',PARAMS{j},',2),[1,3,2]),1);']);
        end

        %interpolate sample to size of B0
        Nx_interB_sample=size(BrainSample.M0,1); Ny_interB_sample=size(BrainSample.M0,2); Nz_interB_sample=size(BrainSample.M0,3);
        x_interB_sample = ([0:Nx_interB_sample-1]-Nx_interB_sample/2+0.5);y_interB_sample = ([0:Ny_interB_sample-1]-Ny_interB_sample/2+0.5);z_interB_sample = ([0:Nz_interB_sample-1]-Nz_interB_sample/2+0.5);
        X_interB_sample=linspace(x_interB_sample(1), x_interB_sample(Nx_interB_sample), matrixsize_signal);
        Y_interB_sample=linspace(y_interB_sample(1), y_interB_sample(Ny_interB_sample), matrixsize_signal);
        Z_interB_sample=linspace(z_interB_sample(1), z_interB_sample(Nz_interB_sample), matrixsize_signal);
        
        for j =1:5
            eval(['BrainSample.',PARAMS{j},'=interp3(y_interB_sample, x_interB_sample, z_interB_sample, BrainSample.',PARAMS{j},', Y_interB_sample'', X_interB_sample, Z_interB_sample, ''spline'');']);
        end
        
        sampleS.M0 = squeeze(BrainSample.M0);
        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(coord_x, 1, []); %x
        sample_straight(3,:) =  reshape(coord_y, 1, []); %y
        sample_straight(4,:) =  reshape(coord_z, 1, []); %z

        %T2
        T2_jemB = BrainSample.T2;
        T2_jemB(T2_jemB<eps) = Inf;    %T2 = 0 is replaced by Inf
        T2_jemB = T2_jemB*10^-3;
        sampleS.T2 = T2_jemB;
        sample_straight(5, :) = reshape(T2_jemB, 1, []); %T2

        %T1
        T1_jemB = BrainSample.T1;
        T1_jemB(T1_jemB<eps) = Inf;    %T1 = 0 is replaced by Inf
        T1_jemB = T1_jemB*10^-3;
        sampleS.T1 = T1_jemB;
        sample_straight(6, :) = reshape(T1_jemB, 1, []); %T1
        
        dB = BrainSample.DB;
        sampleS.DB = dB;
        DB_straight = reshape(squeeze(dB),1,[]);

    case 'Rect1D'
        RectPhantom = rectangularPulse(-settings.reco.FOV/3,settings.reco.FOV/3,x);
          
        sampleS.M0 = squeeze(RectPhantom);
        
        sampleS.T1 = Inf(size(RectPhantom));
        sampleS.T2 = Inf(size(RectPhantom));

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
        dB = zeros(size(RectPhantom));
        sampleS.DB = dB;
        DB_straight = reshape(squeeze(dB),1,[]);
    otherwise
        warning('Not implemented yet!')
end

as(sampleS.M0)
%as(sampleS.M0, 'select', [':,:,', num2str(settings.reco.matrixsize_reco/2)]);
timeElapsed_sample = toc
%clearvars x_reco y_reco x y air_suscept

%% Start Simulation
switch settings.trajectory.type
    case 'ImportFiles'
        [reco_rho_img_func] = simArbFields3D(settings, dB, DB_straight, sampleS, sample_straight);
    case 'Rect1D'
        [reco_rho_img_func] = simArbFields3D(settings, dB, DB_straight, sampleS, sample_straight);
    otherwise
        error('Not implemented!!')
end