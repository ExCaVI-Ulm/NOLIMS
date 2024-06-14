function [settings] = calcMeanB0(settings)
%calcMeanB0 Determine Mean B0-field for all encoding steps -> mixing frequency
%   Input:  -settings struct
%   Output:
%           -settings struct: settings.general.B0 is updated

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

matrixsize_signal_loc = settings.signal.matrixsize_signal;
matrixsize_reco_loc = settings.reco.matrixsize_reco;

Bvec = zeros(settings.trajectory.N_PhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, 4);%1 x, 2 y, 3 z, 4 abs; old: length(time_tot), 
SixRowsLoc = settings.CST.SixRows;

for kkk = 1:settings.trajectory.N_PhaseEnc
    [Bout] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, kkk, matrixsize_signal_loc, matrixsize_reco_loc, settings);
    %Bout(1:3,X,Y,Z); 1:X-Component, 2:Y-Component, 3:Z-Component
    
    Bvec(kkk,:,:,:,1) = squeeze(Bout(1,:,:,:));
    Bvec(kkk,:,:,:,2) = squeeze(Bout(2,:,:,:));
    Bvec(kkk,:,:,:,3) = squeeze(Bout(3,:,:,:));
    Bvec(kkk,:,:,:,4) = sqrt(Bvec(kkk,:,:,:,1).^2 + Bvec(kkk,:,:,:,2).^2 + Bvec(kkk,:,:,:,3).^2);
end

%if only inner part should be used for calculation: meanField_Abs = mean(Bvec(:,(floor(matrixsize_signal_loc/4) : floor(3*matrixsize_signal_loc/4)),(floor(matrixsize_signal_loc/4) : floor(3*matrixsize_signal_loc/4)),(floor(matrixsize_signal_loc/4) : floor(3*matrixsize_signal_loc/4)),4), 'all');
meanField_Abs = mean(Bvec(:,:,:,:,4), 'all');
settings.general.B0 = meanField_Abs;
settings.general.FreqField = 42.577478 * settings.general.B0; %MHz

end

function [Bout] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, jj, matrixsize_signal_loc, matrixsize_reco_loc, settings)
    fileID = fopen([path_rot, filepaths_nonsim{jj}],'r');
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

    BXrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BXrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');
    BYrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BYrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');
    BZrFOV_interpol=interp3(y_Bfield, x_Bfield, z_Bfield, BZrFOV,Y_Bfield', X_Bfield, Z_Bfield, 'linear');
    
    Bout = zeros(3, matrixsize_signal_loc, matrixsize_signal_loc, matrixsize_signal_loc);
    Bout(1,:,:,:) = BXrFOV_interpol;    
    Bout(2,:,:,:) = BYrFOV_interpol;
    Bout(3,:,:,:) = BZrFOV_interpol;
end