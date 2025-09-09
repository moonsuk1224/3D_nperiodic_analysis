%% 3D Reconstruction of Multi-echo FLASH and N-periodic SSFP Data
% This script performs 3D reconstruction of MRI data from raw Siemens k-space (.dat).
% 
% Part 1 - Multi-echo FLASH (reference/ground truth):
%   - coil-by-coil reconstruction
%
% Part 2 - N-periodic SSFP:
%   - Reads 3D phase-cycled SSFP data (N=5 pseudo-states)
%   - Demodulates to obtain F-states (frequency configurations)
%   - Pno POCS
%   - coil-by-coil reconstruction
%   - 
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

%% Save FLASH magnitude images (used as reference)
%save('img_sos_me_flash.mat','img_sos','-v7.3');
%save('img_phase_me_flash.mat','img_phase_flash','-v7.3');
fprintf('FLASH data reconstructed and saved.\n');
fprintf('-----------------------------\n');


%% Export meFLASH Phase to DICOM Format (.IMA)
% DICOM export handles anatomical orientation corrections
cropw=24;
Ny_mid=ceil(Ny/2);
raw_img_crop=zeros(Nx,Ny,Nz,Nc,Necho);
raw_img_crop(:,Ny_mid-cropw:Ny_mid+cropw,:,:,:)=squeeze(d_srt(:,Ny_mid-cropw:Ny_mid+cropw,:,:,:));
raw_img=ifftdim(raw_img_crop,1:3);

for e=1
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

k_fin=zeros(Nx,Ny,Nz,Ncoil,N);
Ny_mid=ceil(Ny/2);
spoiler=20;
cropw=24;
for f=1:N
    fcent=Ny_mid-(f-2-1)*spoiler;
    k_fin(:,Ny_mid-cropw:Ny_mid+cropw-1,:,:,f)=d_fstates_n(:,1+fcent-cropw:fcent+cropw,:,:,f);
end 

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
imagesc(abs(img_phase_sli)); colormap gray; axis image; colorbar;
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

for e=1
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
