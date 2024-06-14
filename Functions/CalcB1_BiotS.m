function [B1_sig, B1] = CalcB1_BiotS(settings,plotFlag)
%CalcB1_BiotS Calculation of B1 of RX array using the Biot-Savart Law
%   Input:  -settings struct
%           -if you wish to plot calculated B1 -> plotFlag = true
%   Output:
%           -B1_sig: calculated B1: matrixsize is that of matrixsize_signal
%           -B1: calculated B1: matrixsize is that of matrixsize_reco


    N = (settings.CoilSens.Ncoil_segments);
    N_c = (settings.CoilSens.NReceiveCoils); 

    Nx_c = settings.signal.matrixsize_signal;
    Ny_c = Nx_c;
    Nz_c = Nx_c;
    
    R_coil = settings.CoilSens.RadiusReceiveCoils;
    R_recArray = settings.CoilSens.RadiusCoilArray;

    curr = settings.CoilSens.Curr;
    u0 = settings.general.u0;
    
    %define points where to evaluate B1
    ax1_px_Coil=linspace(-(settings.reco.FOV)/2, (settings.reco.FOV)/2- settings.trajectory.pixel_width* settings.reco.matrixsize_reco/settings.signal.matrixsize_signal, Nx_c);
    ax2_px_Coil=fliplr(ax1_px_Coil);
    ax3_px_Coil = fliplr(ax1_px_Coil);

    [px_x, px_y, px_z] = meshgrid(ax1_px_Coil,ax2_px_Coil,ax3_px_Coil);
    px_x = permute(px_x, [2 1 3]);
    px_y = permute(px_y, [2 1 3]);
    px_z = permute(px_z, [2 1 3]);
    
    %Describing one coil discretized in N segments: Polar coordinates
    phi_c=(-pi:2*pi/(N-1):pi)'; % angles for one circular coil
    
    %Coil in y-z-plane: Conversion to cartesian coordinates
    X_coil = zeros(N,1) + R_recArray;
    Y_coil = R_coil*cos(phi_c);
    Z_coil = R_coil*sin(phi_c);
    
    B1 = zeros(N_c, Nx_c, Ny_c, Nz_c, 3);

    parfor jj = 1:N_c    
        theta = (jj-1)*(2*pi/N_c);
        %Create RX array by rotating one coil around a rotation axis / here z
        Xcoil_rot = cos(theta)*X_coil - sin(theta)*Y_coil;
        Ycoil_rot = sin(theta)*X_coil + cos(theta)*Y_coil;
        Zcoil_rot = Z_coil; %rotation axis

        L_vec = [Xcoil_rot'; Ycoil_rot'; Zcoil_rot'];  %vector to each segment of the current path

        dl_vec = [diff(L_vec(1,:)); diff(L_vec(2,:)); diff(L_vec(3,:))];  %this vector points towards the next segment -> indicated "direction" of current
        dl_vec = [dl_vec, [Xcoil_rot(1) - Xcoil_rot(end);Ycoil_rot(1) - Ycoil_rot(end);Zcoil_rot(1) - Zcoil_rot(end)]]; %last element must be calculated extra

        for a= 1:Nx_c
            for b=1:Ny_c
                for c = 1:Nz_c

                    r_vec_prime_1pos = [px_x(a,b,c); px_y(a,b,c); px_z(a,b,c)] - L_vec;   % position to evaluate field - position of segment
                    r_mag_prime_1pos = sqrt(sum(r_vec_prime_1pos.^2));    %norm
                    
                    B1_vec = sum(u0*curr/(4*pi)*(cross(dl_vec, r_vec_prime_1pos))./(r_mag_prime_1pos).^3,2);

                    if any(r_mag_prime_1pos<settings.CoilSens.DistMaskWire)      %< settings.CoilSens.DistMaskWire serves as threshold-value to exclude an area around a wire
                        B1_vec = [NaN; NaN; NaN];                                   %exclude wires from calculation
                    end

                    B1x{jj}(a,b,c) =B1_vec(1);
                    B1y{jj}(a,b,c) =B1_vec(2);
                    B1z{jj}(a,b,c) =B1_vec(3);
                end
            end
        end
    end
    
    for j = 1:N_c
        B1(j,:,:,:,1) = B1x{j}(:,:,:);
        B1(j,:,:,:,2) = B1y{j}(:,:,:);
        B1(j,:,:,:,3) = B1z{j}(:,:,:);
    end
    
    clearvars B1_vec r_vec_prime_1pos r_mag_prime_1pos dl_vec L_vec Zc_rot Yc_rot Xc_rot theta Xc Yc Zc phi_c coord_x_pxCoil_dummy px_z px_y px_x ax1_px_Coil ax2_px_Coil ax3_px_Coil
    
    if plotFlag
        as(B1)
    end
    
    %interpolate to matrixsize_reco
    if settings.signal.matrixsize_signal > settings.reco.matrixsize_reco
        Nx_CoilSens=size(B1,2); Ny_CoilSens=size(B1,3); Nz_CoilSens=size(B1,4);
        x_CoilSens = ([0:Nx_CoilSens-1]-Nx_CoilSens/2+0.5);y_CoilSens = ([0:Ny_CoilSens-1]-Ny_CoilSens/2+0.5);z_CoilSens = ([0:Nz_CoilSens-1]-Nz_CoilSens/2+0.5);
        X_CoilSens=linspace(x_CoilSens(1), x_CoilSens(Nx_CoilSens), settings.reco.matrixsize_reco);
        Y_CoilSens=linspace(y_CoilSens(1), y_CoilSens(Ny_CoilSens), settings.reco.matrixsize_reco);
        Z_CoilSens=linspace(z_CoilSens(1), z_CoilSens(Nz_CoilSens), settings.reco.matrixsize_reco);
        B1_sig = B1;
        B1 = zeros(N_c, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, 3);
        for j =1:N_c
            for u = 1:3
                B1(j,:,:,:,u)=interp3(y_CoilSens, x_CoilSens, z_CoilSens, squeeze(B1_sig(j,:,:,:,u)),Y_CoilSens', X_CoilSens, Z_CoilSens, 'linear');
            end
        end
    else
        B1_sig = B1;
    end
end



