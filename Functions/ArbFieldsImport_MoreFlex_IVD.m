function [B_SEM_straight, B_SEM_straight_Reco, BGr_horizontal, BGr_vertikal, BGr_z, BGr_Reco_horizontal, BGr_Reco_vertikal, BGr_Reco_z, B_v_RF, settings] = ArbFieldsImport_MoreFlex_IVD(settings, dB,u)
%ArbFieldsImport_MoreFlex_IVD   Import of encoding field and calculation of the gradients (of the z-component) over each voxel 
%   Input:  -settings struct
%           -3D additional magnetic field caused by susceptibility differences (dB)
%           -u: loop variable
%   Output:
%           -B_SEM_straight: reshaped magnetic field vector used to simulate the signal
%           -B_SEM_straight: reshaped magnetic field vector used for reconstruction
%           -BGr_horizontal: Gradient over each voxel in dimension two for signal IVD
%           -BGr_vertikal: Gradient over each voxel in dimension one for signal IVD
%           -BGr_z: Gradient over each voxel in dimension three for signal IVD
%           -BGr_Reco_horizontal: Gradient over each voxel in dimension two for reconstruction signal model
%           -BGr_Reco_vertikal: Gradient over each voxel in dimension one for signal IVD for reconstruction signal model
%           -BGr_Reco_z: Gradient over each voxel in dimension three for signal IVD for reconstruction signal model
%           -B_v_RF: Magnetic field during RF pulse -> used in Bloch Simulation
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

%IVD
BGr_horizontal = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal^3);
BGr_vertikal = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal^3);
BGr_z = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal^3);

BGr_Reco_horizontal = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco^3);
BGr_Reco_vertikal = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco^3);
BGr_Reco_z = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco^3);

B_SEM_straight = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal^3);
B_SEM_straight_Reco = zeros(1, settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco^3);

B_SEM_z_Comp_Suscept_BlochSim = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);

matrixsize_signal_loc = settings.signal.matrixsize_signal;
matrixsize_reco_loc = settings.reco.matrixsize_reco;

Bvec = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, 4);%1 x, 2 y, 3 z, 4 abs; 
Bvec_Reco = zeros(settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, 4);%1 x, 2 y, 3 z, 4 abs; 
NRAM_loc = settings.general.RAM_StepsPhaseEnc;
SixRowsLoc = settings.CST.SixRows;
FOV_loc = settings.reco.FOV;

%import of files
for kkk = 1:NRAM_loc
    [Bout, Bout_Reco] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, kkk, u, matrixsize_signal_loc, matrixsize_reco_loc, NRAM_loc);
    %Bout(1:3,X,Y,Z); 1:X-Component, 2:Y-Component, 3:Z-Component
    
    Bvec(kkk,:,:,:,1) = squeeze(Bout(1,:,:,:));
    Bvec(kkk,:,:,:,2) = squeeze(Bout(2,:,:,:));
    Bvec(kkk,:,:,:,3) = squeeze(Bout(3,:,:,:));
    Bvec(kkk,:,:,:,4) = sqrt(Bvec(kkk,:,:,:,1).^2 + Bvec(kkk,:,:,:,2).^2 + Bvec(kkk,:,:,:,3).^2);

    Bvec_Reco(kkk,:,:,:,1) = squeeze(Bout_Reco(1,:,:,:));
    Bvec_Reco(kkk,:,:,:,2) = squeeze(Bout_Reco(2,:,:,:));
    Bvec_Reco(kkk,:,:,:,3) = squeeze(Bout_Reco(3,:,:,:));
    Bvec_Reco(kkk,:,:,:,4) = sqrt(Bvec_Reco(kkk,:,:,:,1).^2 + Bvec_Reco(kkk,:,:,:,2).^2 + Bvec_Reco(kkk,:,:,:,3).^2);

    B_SEM_straight(kkk,:) = reshape(squeeze(Bvec(kkk,:,:,:,4)), 1, []);
    B_SEM_straight_Reco(1, kkk, :) = reshape(squeeze(Bvec_Reco(kkk,:,:,:,4)), 1, []);
    B_SEM_z_Comp_Suscept_BlochSim(kkk,:,:,:) = squeeze(Bout(3,:,:,:)) + squeeze((dB));

    %for IVD take only dominant -> z-component, "secular approx"
    Bgr = squeeze(Bvec(kkk,:,:,:,3));
    if settings.general.Suscept
        Bgr = squeeze(B_SEM_z_Comp_Suscept_BlochSim(kkk,:,:,:));
    end

    Bgr_Reco = squeeze(Bvec_Reco(kkk,:,:,:,3));

    %IVD: calculation of gradient over each voxel
    Nx_Bfield_1rot=size(Bgr,1); Ny_Bfield_1rot=size(Bgr,2); Nz_Bfield_1rot=size(Bgr,3);
    x_Bfield_1rot = ([0:Nx_Bfield_1rot-1]-Nx_Bfield_1rot/2+0.5);y_Bfield_1rot = ([0:Ny_Bfield_1rot-1]-Ny_Bfield_1rot/2+0.5);z_Bfield_1rot = ([0:Nz_Bfield_1rot-1]-Nz_Bfield_1rot/2+0.5);
    X_Bfield_1rot=linspace(x_Bfield_1rot(1)-1, x_Bfield_1rot(Nx_Bfield_1rot)+1, Nx_Bfield_1rot+2);
    Y_Bfield_1rot=linspace(y_Bfield_1rot(1)-1, y_Bfield_1rot(Ny_Bfield_1rot)+1, Ny_Bfield_1rot +2);
    Z_Bfield_1rot=linspace(z_Bfield_1rot(1)-1, z_Bfield_1rot(Nz_Bfield_1rot)+1, Nz_Bfield_1rot +2);
	
    Bgr_interp=interp3(z_Bfield_1rot, y_Bfield_1rot, x_Bfield_1rot, squeeze(Bgr),Z_Bfield_1rot', Y_Bfield_1rot', X_Bfield_1rot, 'spline');

    Nx_Bfield_1rot_reco=size(Bgr_Reco,1); Ny_Bfield_1rot_reco=size(Bgr_Reco,2); Nz_Bfield_1rot_reco=size(Bgr_Reco,3);
    x_Bfield_1rot_reco = ([0:Nx_Bfield_1rot_reco-1]-Nx_Bfield_1rot_reco/2+0.5);y_Bfield_1rot_reco = ([0:Ny_Bfield_1rot_reco-1]-Ny_Bfield_1rot_reco/2+0.5);z_Bfield_1rot_reco = ([0:Nz_Bfield_1rot_reco-1]-Nz_Bfield_1rot_reco/2+0.5);
    X_Bfield_1rot_reco=linspace(x_Bfield_1rot_reco(1)-1, x_Bfield_1rot_reco(Nx_Bfield_1rot_reco)+1, Nx_Bfield_1rot_reco+2);
    Y_Bfield_1rot_reco=linspace(y_Bfield_1rot_reco(1)-1, y_Bfield_1rot_reco(Ny_Bfield_1rot_reco)+1, Ny_Bfield_1rot_reco+2);
    Z_Bfield_1rot_reco=linspace(z_Bfield_1rot_reco(1)-1, z_Bfield_1rot_reco(Nz_Bfield_1rot_reco)+1, Nz_Bfield_1rot_reco+2);

    Bgr_interp_Reco=interp3(z_Bfield_1rot_reco, y_Bfield_1rot_reco, x_Bfield_1rot_reco, squeeze(Bgr_Reco),Z_Bfield_1rot_reco',Y_Bfield_1rot_reco', X_Bfield_1rot_reco, 'spline');
    
    Bgr_horz = Bgr_interp(2:size(Bgr_interp,1)-1, :,2:size(Bgr_interp,3)-1);
    Bgr_vert = Bgr_interp(:, 2:size(Bgr_interp,2)-1,2:size(Bgr_interp,3)-1);
    Bgr_z = Bgr_interp(2:size(Bgr_interp,1)-1, 2:size(Bgr_interp,2)-1,:);

    Bgr_Reco_horz = Bgr_interp_Reco(2:size(Bgr_interp_Reco,1)-1, :, 2:size(Bgr_interp_Reco,3)-1);
    Bgr_Reco_vert = Bgr_interp_Reco(:, 2:size(Bgr_interp_Reco,2)-1, 2:size(Bgr_interp_Reco,3)-1);
    Bgr_Reco_z = Bgr_interp_Reco(2:size(Bgr_interp_Reco,1)-1, 2:size(Bgr_interp_Reco,2)-1,:);

    %vertical
    Bgr_vert_m_t = zeros(size(Bgr_vert,1)-2, size(Bgr_vert,2), size(Bgr_vert,3));     		
    Bgr_Reco_vert_m_t = zeros(size(Bgr_Reco_vert,1)-2, size(Bgr_Reco_vert,2), size(Bgr_Reco_vert,3));
    
    for j = 1:size(Bgr_vert,3)
        Bgr_vert_temp_t = squeeze(Bgr_vert(:,:,j));     
        for uu =2:size(Bgr_vert,1)-1
            for kk = 1:size(Bgr_vert,2)
                Bgr_vert_m_t(uu-1,kk,j) = 0.5*(Bgr_vert_temp_t(uu-1,kk) - Bgr_vert_temp_t(uu+1,kk))/(FOV_loc/matrixsize_signal_loc);  %uu-1 switched position with uu+1 since x-axis goes from bottom to top not the other way round: higher value minus lower value coordinate-wise

                if(uu <= size(Bgr_Reco_vert,1)-1 && kk <= size(Bgr_Reco_vert,2) && j <=size(Bgr_Reco_vert,3))
                    Bgr_Reco_vert_temp_t = squeeze(Bgr_Reco_vert(:,:,j));
                    Bgr_Reco_vert_m_t(uu-1,kk,j) = 0.5*(Bgr_Reco_vert_temp_t(uu-1,kk) - Bgr_Reco_vert_temp_t(uu+1,kk))/(FOV_loc/matrixsize_reco_loc);
                end
            end
        end
    end
	
    Bgr_vert_m_diff_enc = Bgr_vert_m_t; 
    Bgr_Reco_vert_m_diff_enc = Bgr_Reco_vert_m_t; 
    		

    %horizontal
    Bgr_horz_m_t = zeros(size(Bgr_horz,1), size(Bgr_horz,2)-2, size(Bgr_horz,3));    
    Bgr_Reco_horz_m_t = zeros(size(Bgr_Reco_horz,1), size(Bgr_Reco_horz,2)-2, size(Bgr_Reco_horz,3));
    
    for j = 1: size(Bgr_horz,3)
        Bgr_horz_temp_t = squeeze(Bgr_horz(:,:,j));                
        for uu =1:size(Bgr_horz,1)
            for kk = 2:size(Bgr_horz,2)-1
                Bgr_horz_m_t(uu,kk-1,j) = 0.5*(Bgr_horz_temp_t(uu,kk+1) - Bgr_horz_temp_t(uu,kk-1))/(FOV_loc/matrixsize_signal_loc);

                if(uu <= size(Bgr_Reco_horz,1) && kk <= size(Bgr_Reco_horz,2)-1&& j <=size(Bgr_Reco_horz,3))
                    Bgr_Reco_horz_temp_t = squeeze(Bgr_Reco_horz(:,:,j));
                    Bgr_Reco_horz_m_t(uu,kk-1,j) = 0.5*(Bgr_Reco_horz_temp_t(uu,kk+1) - Bgr_Reco_horz_temp_t(uu,kk-1))/(FOV_loc/matrixsize_reco_loc);
                end
            end
        end
    end
    Bgr_horz_m_diff_enc = Bgr_horz_m_t;
    Bgr_Reco_horz_m_diff_enc = Bgr_Reco_horz_m_t;
    
        
    %z-dir
    
    Bgr_z_m_t = zeros(size(Bgr_z,1), size(Bgr_z,2), size(Bgr_z,3)-2);    
    Bgr_Reco_z_m_t = zeros(size(Bgr_Reco_z,1), size(Bgr_Reco_z,2), size(Bgr_Reco_z,3)-2);
    
    for j = 1: size(Bgr_z,2)
        Bgr_z_temp_t = squeeze(Bgr_z(:,j,:));
                
        for uu =1:size(Bgr_z,1)
            for kk = 2:size(Bgr_z,3)-1
                Bgr_z_m_t(uu,j, kk-1) = 0.5*(Bgr_z_temp_t(uu,kk+1) - Bgr_z_temp_t(uu,kk-1))/(FOV_loc/matrixsize_signal_loc);

                if(uu <= size(Bgr_Reco_z,1) && kk <= size(Bgr_Reco_z,3)-1&& j <=size(Bgr_Reco_z,2))
                    Bgr_Reco_z_temp_t = squeeze(Bgr_Reco_z(:,j,:));
                    Bgr_Reco_z_m_t(uu,j,kk-1) = 0.5*(Bgr_Reco_z_temp_t(uu,kk+1) - Bgr_Reco_z_temp_t(uu,kk-1))/(FOV_loc/matrixsize_reco_loc);
                end
            end
        end
    end
    Bgr_z_m_diff_enc = Bgr_z_m_t;
    Bgr_Reco_z_m_diff_enc = Bgr_Reco_z_m_t;

    BGr_horizontal(1, kkk,:) = (reshape(Bgr_horz_m_diff_enc, 1, []));
    BGr_vertikal(1, kkk,:) = (reshape(Bgr_vert_m_diff_enc, 1, []));
    BGr_z(1, kkk,:) = (reshape(Bgr_z_m_diff_enc, 1, []));
    
    BGr_Reco_horizontal(1, kkk,:) = (reshape(Bgr_Reco_horz_m_diff_enc, 1, []));
    BGr_Reco_vertikal(1, kkk,:) = (reshape(Bgr_Reco_vert_m_diff_enc, 1, []));
    BGr_Reco_z(1, kkk,:) = (reshape(Bgr_Reco_z_m_diff_enc, 1, []));

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
B_v_RF(:,:,:,:,4) = sqrt(B_v_RF(:,:,:,:,1).^2 + B_v_RF(:,:,:,:,2).^2 + B_v_RF(:,:,:,:,3).^2);

end

function [Bout, Bout_Reco] = ImportDataTXT(filepaths_nonsim, SixRowsLoc, startRow, formatSpec, Nx_CST, Ny_CST, Nz_CST, path_rot, jj,u, matrixsize_signal_loc, matrixsize_reco_loc, NRAM_loc)
    fileID = fopen([path_rot, filepaths_nonsim{jj + (u-1)*NRAM_loc}],'r');
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

    BFOV_reco = sqrt(BXrFOV.^2 + BYrFOV.^2 + BZrFOV.^2);
    BFOV_interpol_Reco=interp3(y_Bfield, x_Bfield, z_Bfield, BFOV_reco,Y_Bfield_Reco', X_Bfield_Reco, Z_Bfield_Reco, 'linear');
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


    Bout = zeros(3, matrixsize_signal_loc, matrixsize_signal_loc, matrixsize_signal_loc);
    Bout(1,:,:,:) = BXrFOV_interpol;
    Bout(2,:,:,:) = BYrFOV_interpol;
    Bout(3,:,:,:) = BZrFOV_interpol;

    Bout_Reco = zeros(3, matrixsize_reco_loc, matrixsize_reco_loc, matrixsize_reco_loc);
    Bout_Reco(1,:,:,:) = BXrFOV_interpol_Reco;
    Bout_Reco(2,:,:,:) = BYrFOV_interpol_Reco;
    Bout_Reco(3,:,:,:) = BZrFOV_interpol_Reco;

end