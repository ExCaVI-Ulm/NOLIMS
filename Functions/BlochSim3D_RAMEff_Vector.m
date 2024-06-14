function [BlochM0] = BlochSim3D_RAMEff_Vector(settings, Sample, Bvec_RF, sample_straight, u)
%BlochSim3D_RAMEff_Vector Bloch Simulation of TX pulse with possibly underlying SEMS taking vector properties into account
%   Input:  -settings struct
%           -sample struct (Sample)
%           -Bvec_RF: Magnetic field during RF pulse
%           -reshaped sample with coordinates, M0, T1, T2 (sample_straight)
%           -u: loop variable
%   Output:
%           -BlochM0: Perpendicular component (wrt B0) of magnetization vector after TX pulse

    gamma_loc = settings.general.gamma;

    spin_r = squeeze(Sample.M0);
    spin_r(spin_r<1E-6) = 0;

    Mlong_init = zeros(settings.general.RAM_StepsPhaseEnc, size(Bvec_RF,2), size(Bvec_RF,3), size(Bvec_RF,4));
    Minit = zeros(settings.general.RAM_StepsPhaseEnc, size(Bvec_RF,2), size(Bvec_RF,3), size(Bvec_RF,4), 3);
    Mvec_pulse = zeros(settings.general.RAM_StepsPhaseEnc,4, size(Bvec_RF,2)*size(Bvec_RF,3)*size(Bvec_RF,4));
    Mvec_pulse_2 = zeros(settings.general.RAM_StepsPhaseEnc,3, size(Bvec_RF,2)*size(Bvec_RF,3)*size(Bvec_RF,4));
    
    dt_Bloch = settings.TX.PulseLength/settings.TX.NumberSamplesBloch;
    time_Bloch = linspace(dt_Bloch, settings.TX.PulseLength, settings.TX.NumberSamplesBloch);
    T2_straight_rot = squeeze(sample_straight(5,:));
    T1_straight_rot = squeeze(sample_straight(6,:));
    
    B_vec_rot = reshape(Bvec_RF, settings.general.RAM_StepsPhaseEnc,settings.signal.matrixsize_signal^3,4); 
    Mean_Babs = zeros(settings.general.RAM_StepsPhaseEnc, 1);
    w = ones(settings.general.RAM_StepsPhaseEnc, length(time_Bloch));  
    
    %Calculation of RF carrier freq + inital thermal magnetization
    for j = 1:settings.general.RAM_StepsPhaseEnc
        Mean_Babs(j) = mean(Bvec_RF(j,:,:,:,4), 'all');
        w(j,:) = settings.general.gamma*Mean_Babs(j)*w(j,:);    %can cause errors if excitation profile should be calculated and displayed: replace Mean_Babs with 3
        Mlong_init(j,:,:,:) = squeeze(Bvec_RF(j,:,:,:,4)).*squeeze(spin_r) / Mean_Babs(j);
        Minit(j,:,:,:,:) = squeeze(Bvec_RF(j,:,:,:,1:3)).*squeeze(spin_r) / Mean_Babs(j);
    end
    Mlong_init_Rot_straight = reshape(Mlong_init, settings.general.RAM_StepsPhaseEnc, size(Bvec_RF,2)*size(Bvec_RF,3)*size(Bvec_RF,4)); 
    Minit_straight = reshape(Minit, settings.general.RAM_StepsPhaseEnc, size(Bvec_RF,2)*size(Bvec_RF,3)*size(Bvec_RF,4), 3);
    %T1 effect
    if settings.general.T1Effect        
        for j = 1:settings.general.RAM_StepsPhaseEnc
            rot_nr = j+(u-1)*settings.general.RAM_StepsPhaseEnc;
            scaling_T1 = settings.sample.T1_scaling;
            for kk = 1:size(Mlong_init_Rot_straight,2)
                Mlong_init_Rot_straight(j,kk) = Mlong_init_Rot_straight(j,kk) * scaling_T1(kk,rot_nr);
                Minit_straight(j,kk,3) = Minit_straight(j,kk,3)* scaling_T1(kk,rot_nr);
            end
        end
    end

    %settings.general.FreqField = settings.general.gamma * mean(Bvec_RF(1,:,:,:,4), 'all')/(2*pi)*10^-6; %in MHz

    %Pulse Definition
    switch settings.TX.pulse_type   %de graaf in vivo spectro
        case 'block'
            B1 = settings.TX.PulseAmpl*ones(1, length(time_Bloch));
        case 'sinc'
            %Five lobe sinc (n=3) 
            B1_norm = sinc(2*3*(time_Bloch - settings.TX.PulseLength/2)/settings.TX.PulseLength);
            B1_max = settings.TX.PulseAmpl;
            B1 = B1_max*B1_norm;
        case 'gaussian10'
            B1_norm = exp(- 0.1 *(2*(time_Bloch - settings.TX.PulseLength/2)/settings.TX.PulseLength)^2);
            B1_max = settings.TX.PulseAmpl;
            B1 = B1_max*B1_norm;
    end
    
   angle_loc =  settings.TX.AngleRF_Theta;

    %B1
    if ~settings.general.B1MapImport
        %if not imported: built now: left circularized uniform B1 lying in
        %x-y plane
        B1x_vec = zeros(settings.general.RAM_StepsPhaseEnc, length(time_Bloch));
        B1y_vec = zeros(settings.general.RAM_StepsPhaseEnc, length(time_Bloch));
        B1z_vec = zeros(settings.general.RAM_StepsPhaseEnc, length(time_Bloch));

        for jj = 1:settings.general.RAM_StepsPhaseEnc
            B1x_vec(jj,1:length(time_Bloch)) = B1.*cos(w(jj,:).*time_Bloch + angle_loc);
            B1y_vec(jj,1:length(time_Bloch)) = -1*B1.*sin(w(jj,:).*time_Bloch + angle_loc);
            B1z_vec(jj,1:length(time_Bloch)) = B1*0;
        end
    B1x_vec = repmat(B1x_vec, [1 1 settings.signal.matrixsize_signal settings.signal.matrixsize_signal settings.signal.matrixsize_signal]);
    B1y_vec = repmat(B1y_vec, [1 1 settings.signal.matrixsize_signal settings.signal.matrixsize_signal settings.signal.matrixsize_signal]);
    B1z_vec = repmat(B1z_vec, [1 1 settings.signal.matrixsize_signal settings.signal.matrixsize_signal settings.signal.matrixsize_signal]);

    else
        B1Import = settings.signal.B1Map;
        B1x_vec = squeeze(B1Import((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc, 1,:,:,:,:));
        B1y_vec = squeeze(B1Import((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc, 2,:,:,:,:));
        B1z_vec = squeeze(B1Import((1:settings.general.RAM_StepsPhaseEnc)+(u-1)*settings.general.RAM_StepsPhaseEnc, 3,:,:,:,:));
    end

    B1x_v_straightSpat = reshape(B1x_vec, settings.general.RAM_StepsPhaseEnc, length(time_Bloch), settings.signal.matrixsize_signal^3); 
    B1y_v_straightSpat = reshape(B1y_vec, settings.general.RAM_StepsPhaseEnc, length(time_Bloch), settings.signal.matrixsize_signal^3); 
    B1z_v_straightSpat = reshape(B1z_vec, settings.general.RAM_StepsPhaseEnc, length(time_Bloch), settings.signal.matrixsize_signal^3); 

   for ll=1:settings.general.RAM_StepsPhaseEnc
            parfor j =1:(size(Bvec_RF,2)*size(Bvec_RF,3)*size(Bvec_RF,4))%par
            
                M_sim = zeros(4,2);                              
                M_sim(:,1) = [squeeze(Minit_straight(ll,j,:)); 1];                                
                          
                for u=2:length(time_Bloch)
                    Btot_z = squeeze(B1z_v_straightSpat(ll,u-1,j)) + squeeze(B_vec_rot(ll,j,3));
                    Btot_y = squeeze(B1y_v_straightSpat(ll,u-1,j)) + squeeze(B_vec_rot(ll,j,2));
                    Btot_x = squeeze(B1x_v_straightSpat(ll,u-1,j)) + squeeze(B_vec_rot(ll,j,1));
                     A = [-1/T2_straight_rot(j),  gamma_loc*Btot_z,  -gamma_loc*Btot_y,       0; ...        
                          -gamma_loc*Btot_z,   -1/T2_straight_rot(j),  gamma_loc*Btot_x ,          0; ...
                           gamma_loc*Btot_y,   -gamma_loc*Btot_x,     -1/T1_straight_rot(j) ,   Mlong_init_Rot_straight(ll,j)/(T1_straight_rot(j)); ...
                            0                     0                            0                        1];
                            
                    M_sim(:,2) = M_sim(:,1) + (A*M_sim(:,1))*dt_Bloch;
                    M_sim(:,1) = M_sim(:,2);                    
                end

                Mvec_pulse = M_sim(:,end); %save M_vector
                %take perpendicular component after Bloch Simulation: only precessing component can be detected
                Mvec_pulse_2(ll,:,j) = squeeze(Mvec_pulse(1:3)) - dot(squeeze(Mvec_pulse(1:3)), squeeze(B_vec_rot(ll,j,1:3))) *squeeze(B_vec_rot(ll,j,1:3)) /norm(squeeze(B_vec_rot(ll,j,1:3)))^2;
            end
   end
    BlochM0 = zeros(settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal,settings.signal.matrixsize_signal,settings.signal.matrixsize_signal,3);
    for jj = 1:3
        BlochM0(:,:,:,:,jj) = reshape(squeeze(Mvec_pulse_2(:,jj,:)), settings.general.RAM_StepsPhaseEnc, settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, []); 
    end
end

