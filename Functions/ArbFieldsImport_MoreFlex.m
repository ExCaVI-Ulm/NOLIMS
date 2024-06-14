function [B_SEM_straight, B_SEM_straight_Reco, B_v_RF, B_v_RF_reco, settings] = ArbFieldsImport_MoreFlex(settings, dB,u)
%ArbFieldsImport_MoreFlex   Import of encoding field
%   Input:  -settings struct
%           -3D additional magnetic field caused by susceptibility differences (dB)
%           -u: loop variable
%   Output:
%           -B_SEM_straight: reshaped magnetic field vector used to simulate the signal
%           -B_SEM_straight: reshaped magnetic field vector used for reconstruction
%           -B_v_RF: Magnetic field during RF pulse -> used in Bloch Simulation
%           -B_v_RF_reco: Magnetic field during RF pulse for reconstruction
%           -settings struct

Nx_CST = settings.CST.Nx_CST;
Ny_CST = settings.CST.Ny_CST;
Nz_CST = settings.CST.Nz_CST;

% search for files to import
path_rot = settings.CST.path_txtFiles;
file_info = natsortfiles(dir(path_rot));
for ooo = 3:(numel(file_info))
    if ~file_info(ooo).isdir    %ignore folders
        filepaths_nonsim{ooo} = file_info(ooo).name;
    end
end
filepaths_nonsim = filepaths_nonsim(~cellfun('isempty',filepaths_nonsim));

for j =1:length(filepaths_nonsim)
    if isempty(regexp(filepaths_nonsim{j},'^d[0-9][0-9]?[0-9]?_', 'once'))
        error('Files do not match required format: "d[degreeOfRot]_"')
    end
end

startRow = 3;

if settings.CST.SixRows
    formatSpec = '%*51s%17f%17f%17f%17f%f%[^\n\r]';
else
    formatSpec = '%*51s%17f%17f%f%[^\n\r]';
end

B_SEM_straight = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal^3);
B_SEM_straight_Reco = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco^3);

B_SEM_z_Comp_Suscept_BlochSim = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);

matrixsize_signal_loc = settings.signal.matrixsize_signal;
matrixsize_reco_loc = settings.reco.matrixsize_reco;

Bvec = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, 4);%1 x, 2 y, 3 z, 4 abs; 
Bvec_Reco = zeros(settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, 4);%1 x, 2 y, 3 z, 4 abs; 
NRAM_loc = settings.general.RAM_StepsPhaseEnc;
SixRowsLoc = settings.CST.SixRows;

for kkk = 1:NRAM_loc
    [Bout, Bout_Reco] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, kkk, u, matrixsize_signal_loc, matrixsize_reco_loc, NRAM_loc, settings);
    %Bout(1:3,X,Y,Z); 1:X-Component, 2:Y-Component, 3:Z-Component
    
    Bvec(kkk,:,:,:,1) = squeeze(Bout(1,:,:,:));
    Bvec(kkk,:,:,:,2) = squeeze(Bout(2,:,:,:));
    Bvec(kkk,:,:,:,3) = squeeze(Bout(3,:,:,:));

    Bvec(kkk,:,:,:,4) = sqrt(Bvec(kkk,:,:,:,1).^2 + Bvec(kkk,:,:,:,2).^2 + Bvec(kkk,:,:,:,3).^2);

    B_SEM_straight(kkk,:) = reshape(squeeze(Bvec(kkk,:,:,:,4)), 1, []);

    Bvec_Reco(kkk,:,:,:,1) = squeeze(Bout_Reco(1,:,:,:));
    Bvec_Reco(kkk,:,:,:,2) = squeeze(Bout_Reco(2,:,:,:));
    Bvec_Reco(kkk,:,:,:,3) = squeeze(Bout_Reco(3,:,:,:));

    Bvec_Reco(kkk,:,:,:,4) = sqrt(Bvec_Reco(kkk,:,:,:,1).^2 + Bvec_Reco(kkk,:,:,:,2).^2 + Bvec_Reco(kkk,:,:,:,3).^2);

    B_SEM_straight_Reco(1, kkk, :) = reshape(squeeze(Bvec_Reco(kkk,:,:,:,4)), 1, []);
    B_SEM_z_Comp_Suscept_BlochSim(kkk,:,:,:) = squeeze(Bout(3,:,:,:)) + squeeze((dB));
end
B_SEM_straight = permute(repmat(B_SEM_straight, 1, 1, 2), [3 1 2]); %keep possibility of time dependence

%magnetic field during Excitation (WITHOUT Exc. magnetic field B1)
B_v_RF = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, 3);
B_v_RF(:,:,:,:,1) = Bvec(:,:,:,:,1);
B_v_RF(:,:,:,:,2) = Bvec(:,:,:,:,2);
B_v_RF(:,:,:,:,3) = Bvec(:,:,:,:,3);
if settings.general.Suscept
    B_v_RF(:,:,:,:,3) = B_SEM_z_Comp_Suscept_BlochSim;
end

B_v_RF_reco = zeros(settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, 3);
B_v_RF_reco(:,:,:,:,1) = Bvec_Reco(:,:,:,:,1);
B_v_RF_reco(:,:,:,:,2) = Bvec_Reco(:,:,:,:,2);
B_v_RF_reco(:,:,:,:,3) = Bvec_Reco(:,:,:,:,3);

B_v_RF(:,:,:,:,4) = sqrt(B_v_RF(:,:,:,:,1).^2 + B_v_RF(:,:,:,:,2).^2 + B_v_RF(:,:,:,:,3).^2);
end

function [Bout, Bout_Reco] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, jj,u, matrixsize_signal_loc, matrixsize_reco_loc, NRAM_loc, settings)
    fileID = fopen(fullfile(path_rot, filepaths_nonsim{jj + (u-1)*NRAM_loc}),'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    fclose(fileID);
    if SixRowsLoc
        BXC1 = dataArray{:, 1};
        BYC1 = dataArray{:, 3};
        BZC1 = dataArray{:, 5};
    else
        BXC1 = dataArray{:, 1};
        BYC1 = dataArray{:, 2};
        BZC1 = dataArray{:, 3};
    end

    f=0; 
    BXrFOV = zeros(Nz_CST, Ny_CST, Nx_CST);
    BYrFOV = zeros(Nz_CST, Ny_CST, Nx_CST);
    BZrFOV = zeros(Nz_CST, Ny_CST, Nx_CST);
    for c=1:Nx_CST
        for b=1:Ny_CST
               for a=1:Nz_CST         
                   f=f+1; 
                   BXrFOV(a,b,c)=BXC1(f);
                   BYrFOV(a,b,c)=BYC1(f);
                   BZrFOV(a,b,c)=BZC1(f);
               end
        end
    end

    %after import CoSy
    %%%%%%%%%%%%%
    % 
    % xC      ________    
    % |     _|______  |
    % v   _|_______ | |
    %    |        | | |
    %    |        | |_|    zC
    %    |        |_|      ^
    %    |________|       /
    %
    % -> yC
    % first dim: x: 1 =neg ->end pos
    %second dim: y: 1= neg -> end pos
    %third dim: z: 1=neg -> end pos
    %%%%%%%%%%%%%

    %!!!!!!!!!!!!
    %field must be in z-dir of the Simulation Coordinate System

    %transform to "new" CoSy
    BXrFOV = flip(BXrFOV,1);    %now x-axis from bottom to top
    BYrFOV = flip(BYrFOV,1);
    BZrFOV = flip(BZrFOV,1);

    %interpolate to desired matrixsize
    
    Nx_Bfield=size(BXrFOV,1); Ny_Bfield=size(BXrFOV,2); Nz_Bfield=size(BXrFOV,3);
    x_Bfield = ([0:Nx_Bfield-1]-Nx_Bfield/2+0.5);y_Bfield = ([0:Ny_Bfield-1]-Ny_Bfield/2+0.5);z_Bfield = ([0:Nz_Bfield-1]-Nz_Bfield/2+0.5);
    X_Bfield=linspace(x_Bfield(1), x_Bfield(Nx_Bfield), matrixsize_signal_loc);
    Y_Bfield=linspace(y_Bfield(1), y_Bfield(Ny_Bfield), matrixsize_signal_loc);
    Z_Bfield=linspace(z_Bfield(1), z_Bfield(Nz_Bfield), matrixsize_signal_loc);

    X_Bfield_Reco=linspace(x_Bfield(1), x_Bfield(Nx_Bfield), matrixsize_reco_loc);
    Y_Bfield_Reco=linspace(y_Bfield(1), y_Bfield(Ny_Bfield), matrixsize_reco_loc);
    Z_Bfield_Reco=linspace(z_Bfield(1), z_Bfield(Nz_Bfield), matrixsize_reco_loc);

    BXrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BXrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');
    BYrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BYrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');
    BZrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BZrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');

    BXrFOV_interpol_Reco=interp3(y_Bfield, x_Bfield, z_Bfield, BXrFOV,Y_Bfield_Reco', X_Bfield_Reco, Z_Bfield_Reco, 'linear');
    BYrFOV_interpol_Reco=interp3(y_Bfield, x_Bfield, z_Bfield, BYrFOV,Y_Bfield_Reco', X_Bfield_Reco, Z_Bfield_Reco, 'linear');
    BZrFOV_interpol_Reco=interp3(y_Bfield, x_Bfield, z_Bfield, BZrFOV,Y_Bfield_Reco', X_Bfield_Reco, Z_Bfield_Reco, 'linear');

    if settings.general.B0MapImport
        BXrFOV_interpol = BXrFOV_interpol + squeeze(settings.signal.B0Map(1,:,:,:));
        BYrFOV_interpol = BYrFOV_interpol + squeeze(settings.signal.B0Map(2,:,:,:));
        BZrFOV_interpol = BZrFOV_interpol + squeeze(settings.signal.B0Map(3,:,:,:));

        if settings.general.B0MapInReco
            BXrFOV_interpol_Reco = BXrFOV_interpol_Reco + squeeze(settings.reco.B0Map(1,:,:,:));
            BYrFOV_interpol_Reco = BYrFOV_interpol_Reco + squeeze(settings.reco.B0Map(2,:,:,:));
            BZrFOV_interpol_Reco = BZrFOV_interpol_Reco + squeeze(settings.reco.B0Map(3,:,:,:));
        end
    end

    %BFOV_reco = sqrt(BXrFOV.^2 + BYrFOV.^2 + BZrFOV.^2);
    %BFOV_interpol_Reco=sqrt(BXrFOV_interpol_Reco.^2 + BYrFOV_interpol_Reco.^2 + BZrFOV_interpol_Reco.^2);%interp3(y_Bfield, x_Bfield, z_Bfield, BFOV_reco,Y_Bfield_Reco', X_Bfield_Reco, Z_Bfield_Reco, 'linear');
    
    Bout = zeros(3, matrixsize_signal_loc, matrixsize_signal_loc, matrixsize_signal_loc);
    Bout(1,:,:,:) = BXrFOV_interpol;    %used for Bloch-Simulations
    Bout(2,:,:,:) = BYrFOV_interpol;
    Bout(3,:,:,:) = BZrFOV_interpol+ settings.general.B0;

    disp('Added Ground field'); %might be omitted if already included in .txt file

    Bout_Reco = zeros(3, matrixsize_reco_loc, matrixsize_reco_loc, matrixsize_reco_loc);
    Bout_Reco(1,:,:,:) = BXrFOV_interpol_Reco;
    Bout_Reco(2,:,:,:) = BYrFOV_interpol_Reco;
    Bout_Reco(3,:,:,:) = BZrFOV_interpol_Reco + settings.general.B0;
end