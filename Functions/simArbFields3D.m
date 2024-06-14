function [reco_rho_img] = simArbFields3D(settings, dB, DB_straight, sampleS, sample_straight)
%simArbFields3D Simulation-File for Arbitrary fields that are imported (3D)
%   Input:  -settings struct
%           -3D additional magnetic field caused by susceptibility differences (dB) & 1D reshaped version (DB_straight)
%           -sample struct (sampleS) & reshaped sample with coordinates, M0, T1, T2 (sample_straight)
%   Output:
%           -reco_rho_img: simulated MR image

    %define time vector
    time_tot = (1:(settings.trajectory.Nsamples))./ settings.trajectory.BW + settings.trajectory.RF_delay;

    time_tot_PE = (1:2*(settings.trajectory.Nsamples))./ settings.trajectory.BW + settings.trajectory.RF_delay;
    
    %initialization of variables needed not only per iteration
    data_meas = zeros(settings.CoilSens.NReceiveCoils, length(time_tot), settings.trajectory.N_PhaseEnc);
    B_SEM_straight_Reco_tot = zeros(1, settings.trajectory.N_PhaseEnc, settings.reco.matrixsize_reco);
    B_SEM_straight_Reco_tot_t = zeros(1, settings.trajectory.N_PhaseEnc, length(time_tot_PE), settings.reco.matrixsize_reco);
    BlochM0_rots = zeros(settings.trajectory.N_PhaseEnc, settings.signal.matrixsize_signal,3);
    CSens_reco_rots = zeros(settings.trajectory.N_PhaseEnc, settings.CoilSens.NReceiveCoils, settings.reco.matrixsize_reco);
    IVDBfield_Reco_tot = zeros(settings.trajectory.N_PhaseEnc, length(time_tot), settings.reco.matrixsize_reco);

    %sample defintion
    rho_vec = repmat(sample_straight(1,:)',[1 settings.trajectory.N_PhaseEnc]);

    %calculate effect of T1 on magnetization that can be used
    if settings.general.T1Effect    
        scaling_T1 = ones(size(rho_vec,1),settings.trajectory.N_PhaseEnc);
        scaling_FA = sin(settings.TX.FlipAngle)*ones(1,settings.trajectory.N_PhaseEnc);
        for j = 2:settings.trajectory.N_PhaseEnc
            for u = 1:size(rho_vec,1)
                scaling_T1(u,j) = 1-(1-squeeze(scaling_T1(u,j-1))*cos(settings.TX.FlipAngle))*exp(-settings.trajectory.TR / sample_straight(6,u));
            end
        end
        rho_vec = scaling_FA.*scaling_T1.*rho_vec;
        settings.sample.T1_scaling = scaling_T1;
    end
    %extent to different coils
    rho_vec = repmat(rho_vec,[1 1 settings.CoilSens.NReceiveCoils]);

    %Coil Sensitivity: calculate B1 once using Biot-Savart, CoilSens is updated for each Encoding field
    if settings.general.CoilSens
        [B1, ~] =CalcB1_BiotS(settings,true);
    else
        %assume uniform B1 in y-direction
        B1 = zeros(settings.CoilSens.NReceiveCoils, settings.signal.matrixsize_signal, 3);
        B1(:,:,2)=1;

        B1_reco = zeros(settings.CoilSens.NReceiveCoils, settings.reco.matrixsize_reco, 3);
        B1_reco(:,:,2)=1;
    end

    tic
    for u = 1:settings.trajectory.N_PhaseEnc / settings.general.RAM_StepsPhaseEnc
        if settings.general.bool_IVD || settings.general.bool_IVD_reco
            %if IVD is simulated by "sincs" ->import fields and calculate the gradient of the z-component of the field over each voxel 
            [B_SEM_straight, B_SEM_straight_Reco, BGr_horizontal, BGr_vertikal, BGr_z, BGr_Reco_horizontal, BGr_Reco_vertikal, BGr_Reco_z, B_v_RF, settings] = ArbFieldsImport_MoreFlex_IVD(settings, dB,u);
        else
            %[B_SEM_straight, B_SEM_straight_Reco, B_v_RF, settings] = ArbFieldsImport_MoreFlex(settings, dB,u);%B_v_RF_reco
            
            %built own encoding fields
            B0_map = zeros(3,settings.signal.matrixsize_signal);

            B0_map(1,:) = 2*sin(20*squeeze(sample_straight(2,:)));
            B0_map(3,:) = settings.general.B0;

            B0_map = settings.general.B0 * B0_map ./(sqrt(settings.general.B0^2 + (2*sin(20*squeeze(sample_straight(2,:)))).^2));
            B_SEM_t = zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot_PE), 3, settings.signal.matrixsize_signal);
            B_SEM_straight = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal);
            B_SEM_straight_t = zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot_PE), settings.signal.matrixsize_signal);
            if u ==1
                B_SEM = B0_map;
                B_SEM_straight(1,:) = settings.general.B0 + settings.trajectory.Gread*squeeze(sample_straight(2,:));
                
                B_SEM1 = B_SEM;
                B_SEM1(3,:) = B_SEM1(3,:) - settings.trajectory.Gread/2*squeeze(sample_straight(2,:)); 
                B_SEM2 = B_SEM;
                B_SEM2(3,:) = B_SEM2(3,:) + settings.trajectory.Gread*squeeze(sample_straight(2,:)); 
                for j = 1:settings.trajectory.Nsamples
                    B_SEM_t(u,j,:,:) = B_SEM1;
                    B_SEM_straight_t(1,j,:) = sqrt((B_SEM1(1,:)).^2 + (B_SEM1(2,:)).^2 + (B_SEM1(3,:)).^2);
                end
                for j = (1+round(settings.trajectory.Nsamples)):size(B_SEM_straight_t,2)
                    B_SEM_t(u,j,:,:) = B_SEM2;
                    B_SEM_straight_t(1,round(j),:) = sqrt((B_SEM2(1,:)).^2 + (B_SEM2(2,:)).^2 + (B_SEM2(3,:)).^2);
                end
            elseif u==2
                B_SEM_straight(1,:) = settings.general.B0 - settings.trajectory.Gread*squeeze(sample_straight(2,:));
            end

            x_reco = fliplr(linspace(-settings.reco.FOV/2, settings.reco.FOV/2- settings.trajectory.pixel_width, settings.reco.matrixsize_reco));

            B_SEM_straight_Reco = zeros(settings.general.RAM_StepsPhaseEnc, settings.reco.matrixsize_reco);
            B_SEM_straight_Reco_t = zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot_PE), settings.reco.matrixsize_reco);
            B_SEM_straight_Reco(1,:) = settings.general.B0 + settings.trajectory.Gread*squeeze(x_reco);

            B0_map_reco = zeros(3,settings.signal.matrixsize_signal);

            B0_map_reco(1,:) = 2*sin(20*squeeze(x_reco));
            B0_map_reco(3,:) = settings.general.B0;

            B0_map_reco = settings.general.B0 * B0_map_reco ./(sqrt(settings.general.B0^2 + (2*sin(20*squeeze(x_reco))).^2));
            B_SEM_reco = B0_map_reco;
            B_SEM1_reco = B_SEM_reco;
            B_SEM1_reco(3,:) = B_SEM1_reco(3,:) - settings.trajectory.Gread/2*squeeze(x_reco); 

            B_SEM2_reco = B_SEM_reco;
            B_SEM2_reco(3,:) = B_SEM2_reco(3,:) + settings.trajectory.Gread*squeeze(x_reco);

            for j = 1:settings.trajectory.Nsamples
                B_SEM_straight_Reco_t(1,j,:) = sqrt((B_SEM1_reco(1,:)).^2 + (B_SEM1_reco(2,:)).^2 + (B_SEM1_reco(3,:)).^2);
            end
            for j = (1+round(settings.trajectory.Nsamples)):size(B_SEM_straight_Reco_t,2)
                B_SEM_straight_Reco_t(1,round(j),:)= sqrt((B_SEM2_reco(1,:)).^2 + (B_SEM2_reco(2,:)).^2 + (B_SEM2_reco(3,:)).^2);
            end

            B_v_RF = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, 4);
            B_v_RF(1,:,1:3) = B0_map.';%here add B0 Map
            B_v_RF(1,:,4) = sqrt(B_v_RF(1,:,1).^2 + B_v_RF(1,:,2).^2 + B_v_RF(1,:,3).^2);

        end

        %Bloch Simulation for TX pulse
        if settings.general.BlochSim
            BlochM0 = BlochSim3D_RAMEff_Vector(settings, sampleS, B_v_RF, sample_straight, u);
            BlochM0_rots((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc,:,:) = BlochM0;

        else
            BlochM0 = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, 3);
            %assume that all initial magnetization lies in y
            BlochM0(:,:,2) = repmat(sampleS.M0, [settings.general.RAM_StepsPhaseEnc, 1]); 
        end

        %which part of magnetization vector can be detected
        M_detect_coil = zeros(settings.CoilSens.NReceiveCoils, settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal);
        for uu = 1:settings.CoilSens.NReceiveCoils
            for ll = 1:settings.general.RAM_StepsPhaseEnc
                for j = 1:settings.signal.matrixsize_signal
                    
                                M_abs = norm(squeeze(BlochM0(ll,j,:)));
                                
                                B1_perp = squeeze(B1(uu,j,:)) - dot(squeeze(B1(uu,j,:)), squeeze(B_SEM2(1:3,j))) * squeeze(B_SEM2(1:3,j)) /(norm(squeeze(B_SEM2(1:3,j)))^2);                                
                                
                                Theta = acos(squeeze(B_v_RF(ll,j,3))/norm(squeeze(B_v_RF(ll,j,1:3))));
                                ez = [0;0;1];
                                if Theta == 0   %B0 in z -> no perpendicular component
                                    Phase_CSens_B1 = atan2(B1(uu,j,2), B1(uu,j,1));
                                    Phase_CSens_M = atan2(BlochM0(ll,j,2), BlochM0(ll,j,1));
                                else
                                    B_perp = squeeze(B_v_RF(ll,j,1:3)) - squeeze(B_v_RF(ll,j,3)) * ez;
                                    rot_axis = cross(ez, B_perp / norm(B_perp));
                                    B1_rot = rot_axis * dot(rot_axis, squeeze(B1(uu,j,:))) + cos(-Theta)*cross(cross(rot_axis, squeeze(B1(uu,j,:))),rot_axis) + sin(-Theta)*cross(rot_axis, squeeze(B1(uu,j,:)));
                                    M_rot = rot_axis * dot(rot_axis, squeeze(BlochM0(ll,j,:))) + cos(-Theta)*cross(cross(rot_axis, squeeze(BlochM0(ll,j,:))),rot_axis) + sin(-Theta)*cross(rot_axis, squeeze(BlochM0(ll,j,:)));
                                    Phase_CSens_B1 = atan2(B1_rot(2), B1_rot(1));
                                    Phase_CSens_M = atan2(M_rot(2), M_rot(1));
                                end
                                
                                if Phase_CSens_B1 < 0
                                    Phase_CSens_B1 = Phase_CSens_B1 + 2*pi;
                                end
                                if Phase_CSens_M < 0
                                    Phase_CSens_M = Phase_CSens_M + 2*pi;
                                end
                                if Phase_CSens_M > Phase_CSens_B1
                                    M_detect_coil(uu,ll,j) = norm(B1_perp)/norm(squeeze(B1(uu,j,:)))*M_abs*exp(1i*acos(dot(squeeze(BlochM0(ll,j,:)), B1_perp) / (norm(B1_perp)*norm(squeeze(BlochM0(ll,j,:))))));
                                else
                                    M_detect_coil(uu,ll,j) = norm(B1_perp)/norm(squeeze(B1(uu,j,:)))*M_abs*exp(-1i*acos(dot(squeeze(BlochM0(ll,j,:)), B1_perp) / (norm(B1_perp)*norm(squeeze(BlochM0(ll,j,:))))));
                                end
                                
                                if isnan(M_detect_coil(uu,ll,j))
                                    M_detect_coil(uu,ll,j) = 0;
                                end
                end
                rho_vec(:,ll + (u-1)*settings.general.RAM_StepsPhaseEnc, uu) = reshape(squeeze(M_detect_coil(uu,ll,:)), 1, []);
            end
        end
        
        %IVD
        %Alternative to Oversampling using matrixsize_signal, matrixsize_reco: simulate IVD per sinc dephasing
        if settings.general.bool_IVD && settings.general.bool_IVD_reco
            IVD = (zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.signal.matrixsize_signal^3));%if used matrixsize_signal = matrixsize_reco
            IVDBfield_Reco = (zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.reco.matrixsize_reco^3)); %no susceptibilities in Fields
            for j = 1:settings.general.RAM_StepsPhaseEnc
                for k = 1:length(time_tot)
                    IVD(j,k,:) = 1/(pi)^3*(abs(pi*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_horizontal(1,j,:))).*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_vertikal(1,j,:)))*pi*pi.*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_z(1,j,:)))));
                    IVDBfield_Reco(j,k,:) = 1/(pi)^3*(abs(pi*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.reco.matrixsize_reco*pi)*time_tot(k)*squeeze(BGr_Reco_horizontal(1,j,:))).*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.reco.matrixsize_reco*pi)*time_tot(k)*squeeze(BGr_Reco_vertikal(1,j,:)))*pi*pi.*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_Reco_z(1,j,:)))));
                end
            end  
            clearvars BGr_horizontal BGr_vertikal BGr_Reco_horizontal BGr_Reco_vertikal BGr_horizontal_Sus_IVDSinc BGr_vertikal_Sus_IVDSinc
        elseif settings.general.bool_IVD
            IVD = (zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.signal.matrixsize_signal^3));%if used matrixsize_signal = matrixsize_reco
            IVDBfield_Reco = (ones(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.reco.matrixsize_reco^3)); %no susceptibilities in Fields
            for j = 1:settings.general.RAM_StepsPhaseEnc
                for k = 1:length(time_tot)
                    IVD(j,k,:) = 1/(pi)^3*(abs(pi*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_horizontal(1,j,:))).*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_vertikal(1,j,:)))*pi*pi.*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_z(1,j,:)))));
                end
            end  
            clearvars BGr_horizontal BGr_vertikal BGr_Reco_horizontal BGr_Reco_vertikal BGr_horizontal_Sus_IVDSinc BGr_vertikal_Sus_IVDSinc
        elseif settings.general.bool_IVD_reco
            IVD = (ones(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.signal.matrixsize_signal^3));%if used matrixsize_signal = matrixsize_reco
            IVDBfield_Reco = (zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.reco.matrixsize_reco^3)); %no susceptibilities in Fields
            for j = 1:settings.general.RAM_StepsPhaseEnc
                for k = 1:length(time_tot)
                    IVDBfield_Reco(j,k,:) = 1/(pi)^3*(abs(pi*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.reco.matrixsize_reco*pi)*time_tot(k)*squeeze(BGr_Reco_horizontal(1,j,:))).*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.reco.matrixsize_reco*pi)*time_tot(k)*squeeze(BGr_Reco_vertikal(1,j,:)))*pi*pi.*sinc(settings.general.gamma*(settings.reco.FOV)/(2*settings.signal.matrixsize_signal*pi)*time_tot(k)*squeeze(BGr_Reco_z(1,j,:)))));
                end
            end  
            clearvars BGr_horizontal BGr_vertikal BGr_Reco_horizontal BGr_Reco_vertikal BGr_horizontal_Sus_IVDSinc BGr_vertikal_Sus_IVDSinc
        else
            IVD = (ones(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.signal.matrixsize_signal));
            IVDBfield_Reco = (ones(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.reco.matrixsize_reco));
        end

        IVDBfield_Reco_tot((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc,:,:) = IVDBfield_Reco;
        
        %% Coil Sensitivity: conversion of B1 to Coil Sensitivity -> depends on SEM -> different for each encoding step
        if settings.general.CoilSens
            [C_rot, C_rot_reco] = CoilSensB1(settings, B1, B_v_RF); 
            C_rot = abs(C_rot);%phase already done with projection onto B1
        else
            C_rot = ones(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils, settings.signal.matrixsize_signal); %uniform Coil Sensitivity
            C_rot_reco = ones(settings.general.RAM_StepsPhaseEnc, settings.CoilSens.NReceiveCoils, settings.reco.matrixsize_reco);
        end
        CSens_reco_rots((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc, :,:) = C_rot_reco;

        %% Larmor freq. in Signal Eq. if low-field
        if settings.general.LowFieldw
            w = settings.general.gamma * squeeze(B_SEM_straight(1,:,:));

            if settings.general.Suscept
                w = w + settings.general.gamma * repmat(DB_straight, [settings.general.RAM_StepsPhaseEnc 1]);
            end
        else
            w = ones(size(B_SEM_straight)); 
        end
        %% Encoding: Calculation of encoding matrix and correponding signal
        phi_enc = zeros(settings.general.RAM_StepsPhaseEnc, length(time_tot), settings.signal.matrixsize_signal);
    
        for jj = 1:settings.general.RAM_StepsPhaseEnc
            
            phi_temp_enc = 1/(2*pi)*settings.general.gamma* cumsum(squeeze(B_SEM_straight_t(jj,:,:) - settings.general.B0), 1)*(time_tot(2)-time_tot(1));
            phi_temp_enc = phi_temp_enc(settings.trajectory.Nsamples+1:end,:);
            if settings.general.T2
                %if T2 should be simulated -> add decaying contribution of transversal magnetization
                phi_temp_enc = phi_temp_enc + 1/(2*pi)*1i*(time_tot./squeeze(sample_straight(5,:)).').';
            end
        
            if settings.general.Suscept
                %if Susceptibility effects should be simulated -> add phase contribution of additional magnet field
                phi_temp_enc = phi_temp_enc + 1/(2*pi)*settings.general.gamma*time_tot.'.*(DB_straight);
            end
            
            phi_enc(jj,:,:) = phi_temp_enc;
        end
    
        clearvars int_IVD_vec phi_temp_enc
        Exp_enc = (exp(2*pi*1i*phi_enc)); %element-wise exponential
        clearvars phi_enc
            
        Nsamples_loc = length(time_tot);
        %calculate signal
        S = (zeros(settings.CoilSens.NReceiveCoils, length(time_tot), settings.general.RAM_StepsPhaseEnc));
        for ll = 1:settings.general.RAM_StepsPhaseEnc
            for mmm=1:settings.CoilSens.NReceiveCoils
                for kk=1:Nsamples_loc
                    S(mmm, kk, ll) = squeeze(w(1,ll,:)).'.*squeeze(C_rot(ll, mmm,:)).'.*squeeze(IVD(ll,kk,:)).'.*squeeze(Exp_enc(ll,kk,:)).'*squeeze(rho_vec(:,ll + (u-1)*settings.general.RAM_StepsPhaseEnc, mmm)); 
                end
                
                if ~settings.general.noise
                    data_meas(mmm, :,(u-1)*settings.general.RAM_StepsPhaseEnc + ll) = squeeze(S(mmm, :,ll));
                else
                    data_meas(mmm, :,(u-1)*settings.general.RAM_StepsPhaseEnc + ll) = (awgn(squeeze(S(mmm, :,ll)), settings.signal.SNR, 'measured'));    %add Noise such that we achieve the desired SNR
                end  
            end
            B_SEM_straight_Reco_tot(:,(u-1)*settings.general.RAM_StepsPhaseEnc + ll,:) = B_SEM_straight_Reco(ll,:);
            B_SEM_straight_Reco_tot_t(:,(u-1)*settings.general.RAM_StepsPhaseEnc + ll,:,:) = B_SEM_straight_Reco_t(ll,:,:);
        end
    end
    signal_gen = toc

    S_cut = data_meas;
    as(BlochM0_rots)
    as(S_cut/8)
    clearvars B_SEM_straight B_SEM_straight_Reco B_v_RF BGr_horizontal BGr_horizontal_Sus_IVDSinc BGr_Reco_horizontal BGr_Reco_vertikal BGr_vertikal BGr_vertikal_Sus_IVDSinc BlochM0 Bvec coord_x coord_x_reco coord_y coord_y_reco coord_z coord_z_reco data_meas dB Exp_enc int_IVD_vec phi_enc S vec vec_reco rho_vec sampleS
    %% Reconstruction
    disp('Reco');tic;
    
    %if CoilSens is used -> transform S_cut
    if settings.general.CoilSens
        Signal_concat_coil = zeros(settings.trajectory.N_PhaseEnc*settings.trajectory.Nsamples, settings.CoilSens.NReceiveCoils);
        for u = 1:settings.CoilSens.NReceiveCoils            
            signal_coil_tmp = squeeze(S_cut(u,:,:));
            Signal_concat_coil(:,u) = signal_coil_tmp(:);
        end
        S_cut = Signal_concat_coil;
    end

    amplitude_mod1 = zeros(settings.CoilSens.NReceiveCoils, settings.trajectory.N_PhaseEnc, settings.reco.matrixsize_reco);
    %amplitude modulation due to sensitivity differences
    for uu = 1:settings.CoilSens.NReceiveCoils
        for ll = 1:settings.trajectory.N_PhaseEnc
            for j = 1:settings.reco.matrixsize_reco
                B1_perp_r = squeeze(B1_reco(uu,j,:)) - dot(squeeze(B1_reco(uu,j,:)), squeeze(B_SEM2_reco(1:3,j))) * squeeze(B_SEM2_reco(1:3,j)) /(norm(squeeze(B_SEM2_reco(1:3,j)))^2);                                
                amplitude_mod1(uu,ll,j) = norm(B1_perp_r)/norm(squeeze(B1_reco(uu,j,:)));
            end
        end
    end
    amplitude_mod1((amplitude_mod1 == 0)) = 1;    %if 0 -> no sensitivity -> not possible to recover anything, also makes problems with inversion of E

    %for Reco: setting up the encoding matrix either for PCG or ART (RAMSavingReco)
    if ~settings.general.RAMSavingReco
        for j =1:settings.trajectory.N_PhaseEnc
                        
            phi_temp = 1/(2*pi)*settings.general.gamma* cumsum(squeeze(B_SEM_straight_Reco_tot_t(1,j,:,:) - settings.general.B0), 1)*(time_tot(2)-time_tot(1));
            phi_temp = phi_temp(settings.trajectory.Nsamples+1:end,:);

            phi_vert{j} = phi_temp;

            if settings.general.bool_IVD_reco
                %if IVD is incorporated into the Reconstruction Model
                IVDBfield_Reco_tot_vert{j} = squeeze(IVDBfield_Reco_tot(j,:,:));
            end

            if settings.general.LowFieldw
                %if LowField -> w(r) is taken into account in the signal equation
                w_Reco_time{j} = settings.general.gamma * repmat(squeeze(B_SEM_straight_Reco_tot(1,j,:)).', [settings.trajectory.Nsamples, 1]);
            end
        end
        phi_reco_Nrot = vertcat(phi_vert{:});
        clearvars phi_vert phi_temp
        E_reco = (exp(1i*2*pi*phi_reco_Nrot));
        
        amplitude_mod = repmat(squeeze(amplitude_mod1), [1, round(settings.trajectory.Nsamples)]).';
        amplitude_mod((amplitude_mod == 0)) = 1;    %if 0 -> no sensitivity -> not possible to recover anything, also makes problems with inversion of E
        E_reco = amplitude_mod .*E_reco;

        if settings.general.bool_IVD_reco
            IVDBfield_Reco_tot_Nrot = vertcat(IVDBfield_Reco_tot_vert{:});
            E_reco = IVDBfield_Reco_tot_Nrot.*E_reco;
            clearvars IVDBfield_Reco_tot IVDBfield_Reco_tot_vert IVDBfield_Reco_tot_Nrot
        end  

        if settings.general.LowFieldw
            w_Reco = vertcat(w_Reco_time{:});
            E_reco = w_Reco.*E_reco;
            clearvars w_Reco w_Reco_time
        end

        if settings.general.CoilSens
            %structure of E_reco
% 
%             coil loop
%                 encoding loop
%                     time loop
% 
%                     end time loop
%                 end encoding loop
%             end coil loop

            E_reco = repmat(E_reco,[settings.CoilSens.NReceiveCoils, 1]);
            for u = 1:settings.CoilSens.NReceiveCoils
                
                %for each rotation: add time steps and concate them
                for kk = 1:settings.trajectory.N_PhaseEnc
                    tmp_timeConcat{kk} = repmat(squeeze(CSens_reco_rots(kk,u,:)).', [settings.trajectory.Nsamples, 1]);
                end
                timeConcat_full = vertcat(tmp_timeConcat{:});

                %then concat for all RX coils
                if u == 1
                    Csens_Concat = timeConcat_full;
                else
                    tmp_Csens_Concat = timeConcat_full;
                    Csens_Concat = [Csens_Concat; tmp_Csens_Concat];
                end     

            end
            E_reco = Csens_Concat.*E_reco;
        end

        
        iteration=100;
    
        S_backStraightCut = S_cut(:);
        reco_rho_straight = pcg(E_reco'*E_reco, E_reco'*S_backStraightCut, 1e-10, iteration, diag(diag(E_reco'*E_reco)));
        reco_rho_img = reco_rho_straight;
    else
        %Kaczmarz method /ART reco -> calculate each row of E individually -> reduction of RAM demands
        
        %Parameter that worked well
        max_it_ART = 360;
        lambda = 0.08; 
        N_ART = settings.reco.matrixsize_reco;
        x0_ART = zeros(N_ART,1);
        S_backStraightCut = S_cut(:);

        %if your GPU VRAM is large enough -> do the reconstruction on the GPU
        if settings.general.gpuComp
            B_SEM_straight_Reco_tot = gpuArray((B_SEM_straight_Reco_tot_t));
            time_tot = gpuArray(time_tot);
            S_backStraightCut = gpuArray(S_backStraightCut);
            X = gpuArray(x0_ART);
            IVDBfield_Reco_tot = gpuArray(IVDBfield_Reco_tot);
            CSens_reco_rots = gpuArray(CSens_reco_rots);
            amplitude_mod1 = gpuArray(amplitude_mod1);
        else
            X = (x0_ART);
        end

        clearvars x0_ART
        for j = 1:max_it_ART
            for jj=1:length(S_backStraightCut)

                %calculation of current coil_number / encoding step / time index from jj
                if ~mod(jj, length(time_tot)*settings.trajectory.N_PhaseEnc)
                   c_index =  floor(jj/(length(time_tot)*settings.trajectory.N_PhaseEnc));
                else
                    c_index = floor(jj/(length(time_tot)*settings.trajectory.N_PhaseEnc))+1;
                end
                
                if ~mod(jj,length(time_tot))
                    t_index =length(time_tot);
                else
                    t_index = mod(mod(jj, length(time_tot)*settings.trajectory.N_PhaseEnc), length(time_tot));
                end
                
                rot_index = floor((jj - (c_index-1) * length(time_tot)*settings.trajectory.N_PhaseEnc - 1)/length(time_tot))+1;
                int_B = cumsum(squeeze(B_SEM_straight_Reco_tot(1,rot_index,:,:) - settings.general.B0), 1);
                int_B = int_B(round(settings.trajectory.Nsamples)+1:end,:);
                An = exp(1i*settings.general.gamma*int_B(t_index,:)*(time_tot(2)-time_tot(1)));                

                if settings.general.bool_IVD_reco
                    An = squeeze(IVDBfield_Reco_tot(rot_index,t_index,:)).'.*An;
                end

                if settings.general.CoilSens
                    An = squeeze(CSens_reco_rots(rot_index, c_index,:)).'.*An;
                end

                if settings.general.LowFieldw
                    An = settings.general.gamma*squeeze(B_SEM_straight_Reco_tot(1,rot_index,t_index,:)).'.*An;
                end

                %amplitude-correction
                An = squeeze(amplitude_mod1(c_index,rot_index,:)).'.*An;

                unitLen = norm(An);
                d = (S_backStraightCut(jj)-An*X)./unitLen^2;
                X = X+(lambda*d.*conj(An)).';
            end
        end
        
        if settings.general.gpuComp
            reco_rho_img = gather(X);
        else
            reco_rho_img = X;
        end

    end

    as(reco_rho_img, 'Img Func');
    elapsedtime_reco = toc
end