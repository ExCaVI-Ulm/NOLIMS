function [C_rot_resh, C_rot_resh_reco] = CoilSensB1(settings,B1, Bvec)
%CoilSensB1 Calculate Coil Sensitivity of RX Array with respect to B1 calculated with Biot-Savart
%   Input:  -settings struct
%           -B1: Magnetic field produced by RX array calculated with the Biot-Savart Law
%           -Bvec: Encoding field / B0
%   Output:
%           -C_rot_resh: Reshaped Signal Coil Sensitivity for encodings "1:RAM_loc"
%           -C_rot_reco: Reco Coil Sensitivity for encodings "1:RAM_loc", interpolated from C_rot

%Magnitude
Mag_CSens = zeros(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);

%Phase
Phase_CSens = zeros(size(Mag_CSens));

%CoilSens
C_rot = zeros(size(Mag_CSens));
C_rot_resh = zeros(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils, settings.signal.matrixsize_signal^3);
C_rot_resh_reco = zeros(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils, settings.reco.matrixsize_reco^3);

Nx_CoilSens=size(C_rot,3); Ny_CoilSens=size(C_rot,4); Nz_CoilSens=size(C_rot,5);
x_CoilSens = ([0:Nx_CoilSens-1]-Nx_CoilSens/2+0.5);y_CoilSens = ([0:Ny_CoilSens-1]-Ny_CoilSens/2+0.5);z_CoilSens = ([0:Nz_CoilSens-1]-Nz_CoilSens/2+0.5);
X_CoilSens=linspace(x_CoilSens(1), x_CoilSens(Nx_CoilSens), settings.reco.matrixsize_reco);
Y_CoilSens=linspace(y_CoilSens(1), y_CoilSens(Ny_CoilSens), settings.reco.matrixsize_reco);
Z_CoilSens=linspace(z_CoilSens(1), z_CoilSens(Nz_CoilSens), settings.reco.matrixsize_reco);
tic
max_coilsens = zeros(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils);

for m = 1:settings.general.RAM_StepsPhaseEnc
    for ll = 1:settings.CoilSens.NReceiveCoils
        for j = 1:settings.signal.matrixsize_signal
            for u = 1: settings.signal.matrixsize_signal
                for k = 1:settings.signal.matrixsize_signal
                    Mag_CSens(m, ll, j,u,k) = norm(squeeze(B1(ll,j,u,k,:)) - dot(squeeze(B1(ll,j,u,k,:)), squeeze(Bvec(m,j,u,k,1:3)))/norm(squeeze(Bvec(m,j,u,k,1:3)))^2 * squeeze(Bvec(m,j,u,k,1:3)));
                    
                    Theta = acos(squeeze(Bvec(m,j,u,k,3))/norm(squeeze(Bvec(m,j,u,k,1:3))));
                    ez = [0;0;1];
                    if Theta == 0   %B0 in z -> no perpendicular component
                        Phase_CSens(m, ll, j,u,k) = atan2(-B1(ll,j,u,k,2), B1(ll,j,u,k,1));
                    else
                        B_perp = squeeze(Bvec(m,j,u,k,1:3)) - squeeze(Bvec(m,j,u,k,3)) * ez;
                        rot_axis = cross(ez, B_perp / norm(B_perp));
                        B1_rot = rot_axis * dot(rot_axis, squeeze(B1(ll,j,u,k,:))) + cos(-Theta)*cross(cross(rot_axis, squeeze(B1(ll,j,u,k,:))),rot_axis) + sin(-Theta)*cross(rot_axis, squeeze(B1(ll,j,u,k,:)));
                        Phase_CSens(m, ll, j,u,k) = atan2(-B1_rot(2), B1_rot(1));
                    end
                    C_rot(m,ll,j,u,k) = Mag_CSens(m, ll, j,u,k)*exp(1i*Phase_CSens(m, ll, j,u,k));
                    
                    if(Mag_CSens(m, ll, j,u,k) > max_coilsens(m,ll))
                        max_coilsens(m,ll) = Mag_CSens(m, ll, j,u,k);
                    end
                end
            end
        end
        C_rot(m,ll,:,:,:) = C_rot(m,ll,:,:,:) / max_coilsens(m,ll);

        %set coilsens to one where B1 was NaN -> wire
        tmp_coilsens = squeeze(C_rot(m,ll,:,:,:));
        tmp_coilsens(isnan(tmp_coilsens)) = 1;
        C_rot(m,ll,:,:,:) =tmp_coilsens;

        C_rot_resh(m,ll,:) = reshape(squeeze(C_rot(m,ll,:,:,:)), 1, []);
        %Reco C_Sens
        C_rot_resh_reco(m,ll,:) = reshape(interp3(y_CoilSens, x_CoilSens, z_CoilSens, squeeze(C_rot(m,ll,:,:,:)),Y_CoilSens', X_CoilSens, Z_CoilSens, 'linear'), 1, []);
    end
end
toc
as(C_rot)
as(C_rot_resh_reco)


end



