%% 3D Reconstruction of Multi-echo FLASH and N-periodic SSFP Data
% This script performs 3D reconstruction of MRI data from raw Siemens k-space (.dat).
% 
% Part 1 - Multi-echo FLASH (reference/ground truth):
%   - coil-by-coil reconstruction
%
% Part 2 - N-periodic SSFP:
%   - Reads 3D phase-cycled SSFP data (N=5 pseudo-states)
%   - Demodulates to obtain F-states (frequency configurations)
%   - Performs 2D POCS partial Fourier reconstruction per slice and coil
%   - coil-by-coil reconstruction
%   - write DICOM
% Note:
% 

clear; clc; close all;
parpool;

%% Add external functions to path and navigate to data folder
startdir=pwd;
% Update these paths if they are different on your system
addpath(genpath("..\external"))
cd('../data')

%% ===== MULTI-ECHO FLASH RECONSTRUCTION (3D ACQUISITION) =====
fprintf('--- Processing 3D Multi-echo FLASH Data ---\n');

% Read raw data file 254 (3D FLASH sequence)
twix_flash = mapVBVD(254);
twix_flash.image.flagRemoveOS = true;  % Remove 2x oversampling in readout direction

% Get raw k-space data
d_orig = twix_flash.image();

% Reorder dimensions from Siemens format to standard ordering
% Original: [Col, Cha, Lin, Par, ..., Eco] 
% Target:   [RO, PE, SL, CH, ECHO] for [128, 128, 128, 32, 5]
d_srt = permute(d_orig, [1 3 4 2 8 5 6 7]);
d_srt = squeeze(d_srt);
fprintf('FLASH data loaded: 128x128x128, 32 coils, 5 echoes\n');

% Extract dimensions
[Nx, Ny, Nz, Nc, Necho] = size(d_srt);

% % Reconstruct each echo using 2D FFT per slice 
% img_mag = zeros(Nx, Ny, Nz, Necho);
% for echo = 1:Necho
%     for slc = 1:Nz
%         % Apply 2D inverse FFT to each slice, then combine coils
%         img_sos(:,:,slc,echo) = sqrt(sum(abs(ifft2c(squeeze(d_srt(:,:,slc,:,echo)))).^2,3));
%     end
% end


%% Save FLASH magnitude images (used as reference)
%save('img_sos_me_flash.mat','img_sos','-v7.3');
%save('img_phase_me_flash.mat','img_phase_flash','-v7.3');
fprintf('FLASH data reconstructed and saved.\n');
fprintf('-----------------------------\n');


%% Export meFLASH Phase to DICOM Format (.IMA)
% DICOM export handles anatomical orientation corrections
raw_img=ifftdim(d_srt,1:3);
for e=1:Nc
    img_phase_flash = squeeze(angle(raw_img(:,:,:,e,:)));
    img_mag=squeeze(abs(raw_img(:,:,:,e,:)));
    fname = sprintf('img_phase_meFLASH_coil%02d.mat', e);
    save(fname,'img_phase_flash','-v7.3');
    fname = sprintf('img_mag_meFLASH_coil%02d.mat', e);
    save(fname,'img_mag','-v7.3');
    
    TE1_ms = 4;                   % Echo time at center of readout
    voxelSize_mm = [2.0 2.0 2.0]; % Isotropic 2mm voxels
    B0_T = 3;                      % 3 Tesla scanner
    TR_ms = 41;
    deltaTE_ms = 8;
    patientName = 'Volunteer^FLASH_3D';
    patientID = 'FLASH001';
    
    % Create output directory structure
    prefix = sprintf('out/tr%02dms_coil%d', TR_ms, e);
    if ~exist('out', 'dir')
        mkdir('out');
    end

    % Export with orientation correction (handled in export function)
    % Creates separate DICOM series for each echoes
    meflash2dcm(prefix, img_mag, img_phase_flash, TE1_ms, deltaTE_ms, TR_ms, voxelSize_mm, B0_T, patientName, patientID, 001);
end

%% ===== N-PERIODIC SSFP RECONSTRUCTION (3D PARTIAL FOURIER) =====
% Acquisition parameters for 3D phase-cycled SSFP:
%   - N = 5 pseudo-states (phase cycling steps)
%   - TR = 8.00 ms, TE = 4.00 ms, Flip angle = 8°
%   - Central configuration: F2 (best SNR)
%   - Gradient spoiler shift: 20 PE lines between states
%   - 3D matrix ~ 128x128x128
%   - Partial Fourier in PE (recovered with POCS)

file_name = 255;  % Raw data file number for SSFP
N = 5;            % Number of phase-cycled pseudo-states
TR = 8;           % Repetition time in ms
cfg = 2;          % Central F-state configuration (F2 has symmetric sampling)
spl_lines = 20;   % Gradient spoiler induced PE line shift per F-state

fprintf('--- Processing 3D n-periodic SSFP (N=%d, cfg=F%d) ---\n', N, cfg);

% Read SSFP raw data
twix = mapVBVD(file_name);
twix.image.flagRemoveOS = true;  % Remove readout oversampling
d_orig = squeeze(twix.image());
d_pss_n = permute(d_orig, [1 3 4 2 5]);  % Final shape: [128,128,128,32,5]
[Nx,Ny,Nz,Ncoil,N] = size(d_pss_n);
fprintf('SSFP dataset: %dx%dx%d volume, %d coils, %d states\n', Nx, Ny, Nz, Ncoil, N);

%% Demodulate Phase-Cycled Data to F-states
% Each F-state represents a different frequency configuration
% Demodulation: multiply by exp(-i*2π*f*n/N) where f=state index, n=cycle index
fprintf('Demodulating to F-states...\n');
fstates = 0:N-1;
d_fstates_n = zeros(Nx,Ny,Nz,Ncoil,N, 'single');
i=1;
for f_idx = fstates
    % Apply phase demodulation and average across cycles
    demod_factor = exp(f_idx*1i*2*pi*(0:N-1)/N);
    demod_factor = reshape(demod_factor, 1, 1, 1, 1, []);
    d_fstates_n(:,:,:,:,i) = mean(d_pss_n .* demod_factor, 5);
    i=i+1;
end

%% Quick visualisation
figure(1)
img_disp=ifftdim(d_fstates_n,1:3);
sli=ceil(Nz/2);
img_disp_sli=squeeze(img_disp(:,:,sli,:,:));
img_disp_sos = squeeze(sqrt(sum(abs(img_disp_sli).^2,3)));
img_disp_sos=reshape(img_disp_sos,Nx,[]);
%img_disp_sos=rot90(img_disp_sos);
imagesc(img_disp_sos); colormap gray; axis image;
set(gca, 'XTick', [], 'YTick', []);
xlabel('From left to right: F0, F1, F2, F3, F4, respectively');
title('Initial reconstruction (before POCS)');

%% 2D POCS Partial Fourier Reconstruction (slice-by-slice)
% POCS = Projection Onto Convex Sets
% Used to recover missing k-space data from partial sampling
% Each F-state has different k-space coverage pattern due to gradient spoiling

iters=10; % choose number of POCS iterations
k_fin=zeros(Nx,Ny,Nz,Ncoil,N);
    %Loop through each of the N F-states
for f_idx=1:N
    for sli=1:Nz
        S = zeros(Nx,Ny,Ncoil); % always 2D to align with Pete's code
        %Work out which lines we sample, and which we don't (which POCS will
        %calculate)
        pe_mask = circshift([ones(1,Ny),zeros(1,Ny)],spl_lines*(f_idx-cfg-1));
        pe_lines = find(pe_mask(1:Ny));
        non_pe_lines = setxor(pe_lines,1:Ny);

        %Choose the central lines around each echo to work out the image phase
        mid_mask = zeros(size(S,[1 2]));
        sym_sz = sum(pe_mask(1:Ny))-(Ny/2)-1;
        mid_mask(:,(Ny/2)-sym_sz:(Ny/2)+sym_sz-1)=1;
        n_lines = length(pe_lines);

        %Fill in measured lines into S
        if f_idx>(cfg+1)  % Is F-state to left or right of centre?
            S(:,pe_lines,:)=squeeze(d_fstates_n(:,1:n_lines,sli,:,f_idx));    
        else 
            S(:,pe_lines,:)=squeeze(d_fstates_n(:,end-n_lines+1:end,sli,:,f_idx));    
        end

        % Compute phase consistency term for POCS
        phase_im = exp(1i.*angle(ifftdim(mid_mask.*S,1:2)));

        %% Run POCS reconstruction for each coil and visualise results
        for idx_c=1:Ncoil

            % fill k-space with known data
            k_pocs = zeros(size(S,[1 2]));             
            k_pocs(:,pe_lines)=S(:,pe_lines,idx_c); 

            for it = 1:iters   %for each POCS iteration

                im_pocs = abs(ifftdim(k_pocs,1:2)).*phase_im(:,:,idx_c);    % enforce phase consistency
                k_new = fftdim(im_pocs,1:2);                                % update k-space guess
                k_pocs(:,non_pe_lines)=k_new(:,non_pe_lines);               % combine with known data

            end %ending iteration loop for POCS

            % store final k-space results per coil
            k_fin(:,:,sli,idx_c,f_idx)=single(k_pocs);
    
        end  %ending coil loop for POCS
    end
end %ending F-state loop for this dataset

%% Sanity check for POCS (in k-space)

figure(2)
img_pocs=ifftdim(squeeze(k_fin),1:3);
sli=ceil(Nz/2);
img_pocs_sos = squeeze(sqrt(sum(abs(img_pocs).^2,4)));
c=1;
img_pocs_sli=squeeze(img_pocs_sos(:,:,sli,:));
img_pocs_sli=reshape(img_pocs_sli,Nx,[]);
img_phase_sli=squeeze(img_pocs(:,:,sli,c,:));
img_phase_sli=reshape(img_phase_sli,Nx,[]);
imagesc(angle(img_phase_sli)); colormap jet; axis image; colorbar;
set(gca, 'XTick', [], 'YTick', []);
xlabel('From left to right: F0, F1, F2, F3, F4, respectively');
title('Sanity Check for POCS (in k-space)');



%% Extract Magnitude for Each F-state
fprintf('Creating Mag images across F-states...\n');
img_mag_ssfp = zeros(Nx, Ny, Nz, N, 'single');    % Magnitude images
img_mag_ssfp=abs(img_pocs);

%% Extract Phase for Each F-state
fprintf('Creating Phase images across F-states...\n');
img_phase_ssfp = zeros(Nx, Ny, Nz, N, 'single');    % Magnitude images
img_phase_ssfp=angle(img_pocs);

%% Export meFLASH Phase to DICOM Format (.IMA)
% DICOM export handles anatomical orientation corrections

for e=1:Nc
    img_phase_coil = squeeze(img_phase_ssfp(:,:,:,e,:));
    img_mag_coil=squeeze(img_mag_ssfp(:,:,:,e,:));
    fname = sprintf('img_phase_ssfp_coil%02d.mat', e);
    save(fname,'img_phase_coil','-v7.3');
    fname = sprintf('img_mag_ssfp_coil%02d.mat', e);
    save(fname,'img_mag_coil','-v7.3');
    
    TE0_ms = 4;                   % Echo time at center of readout
    voxelSize_mm = [2.0 2.0 2.0]; % Isotropic 2mm voxels
    B0_T = 3;                      % 3 Tesla scanner
    patientName = 'Volunteer^SSFP_3D';
    patientID = 'SSFP001';
    
    % Create output directory structure
    prefix = sprintf('out/tr%02dms_coil%d', TR, e);
    if ~exist('out', 'dir')
        mkdir('out');
    end
    
    % Export with orientation correction (handled in export function)
    % Creates separate DICOM series for each F-state
    ssfp2dcm(prefix, img_mag_coil, img_phase_coil, TE0_ms, TR, voxelSize_mm, B0_T, patientName, patientID, 001);
end

%% Save MATLAB Files

fprintf('3D SSFP reconstruction complete (N=%d F-states).\n', N);
fprintf('All datasets processed with correct orientation.\n');

% Clean up parallel pool
delete(gcp('nocreate'));
