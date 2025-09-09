function meflash2dcm(prefix, img_mag, img_phase, TE1_ms, deltaTE_ms, TR_ms, voxelSize_mm, B0_T, patientName, patientID, seriesNumberBase)
% Export multi-echo FLASH magnitude and phase images to DICOM with corrected orientation
% Siemens-like phase scaling (0..4095 unsigned) with slope/intercept -> [-pi, +pi)
%
% prefix         : output folder prefix (created if missing), e.g. 'FLASH_out/tr41ms_e5'
% img_mag/phase  : single, size [Nx,Ny,Nz,Necho]; phase in radians (-pi..pi)
% TE1_ms         : first echo time in ms
% deltaTE_ms     : echo spacing in ms
% TR_ms          : repetition time in ms
% voxelSize_mm   : [dx dy dz] in mm
% B0_T           : Tesla
% patientName/ID : strings
% seriesNumberBase : int base for SeriesNumber (mag uses base, phase uses base+1)

arguments
    prefix (1,:) char
    img_mag single
    img_phase single
    TE1_ms (1,1) double = 4
    deltaTE_ms (1,1) double = 8
    TR_ms (1,1) double = 41
    voxelSize_mm (1,3) double = [2.0 2.0 2.0]
    B0_T (1,1) double = 3
    patientName (1,:) char = 'Anon^meFLASH'
    patientID   (1,:) char = '000000'
    seriesNumberBase (1,1) double = 0
end

% ---------- shapes ----------
sz = size(img_mag); if numel(sz)==3, sz(4)=1; end
[Nx,Ny,Nz,Necho] = deal(sz(1),sz(2),sz(3),sz(4));
assert(isequal(size(img_phase), [Nx,Ny,Nz,Necho]), 'mag/phase shape mismatch');

% Echo times
TE_ms = TE1_ms + (0:Necho-1) * deltaTE_ms;
fprintf('meFLASH with %d echoes: TE = [', Necho); fprintf(' %.1f', TE_ms); fprintf(' ] ms\n');

%% ---------- orientation correction (simple flips to match typical Siemens axial) ----------
fprintf('Applying orientation corrections...\n');
img_mag_c   = zeros(size(img_mag),   'single');
img_phase_c = zeros(size(img_phase), 'single');
for e = 1:Necho
    tmpm = flip(flip(img_mag(:,:,:,e), 2), 3);   % flip L-R (dim2) and slice (dim3)
    tmpp = flip(flip(img_phase(:,:,:,e), 2), 3);
    img_mag_c(:,:,:,e)   = tmpm;
    img_phase_c(:,:,:,e) = tmpp;
end
img_mag   = img_mag_c;
img_phase = img_phase_c;

% ---------- make output dirs ----------
outMagDir = fullfile(prefix,'DICOM_MAG');   if ~exist(outMagDir,'dir'),   mkdir(outMagDir);   end
outPhsDir = fullfile(prefix,'DICOM_PHASE'); if ~exist(outPhsDir,'dir'), mkdir(outPhsDir); end

% ---------- scaling ----------
% Magnitude: uint16 with linear rescale back to native intensity
mag_max   = max(img_mag(:));  if mag_max<=0, mag_max = 1; end
mag_scale = 65535 / double(mag_max);
mag_rescale_slope     = 1 / mag_scale;
mag_rescale_intercept = 0;

% Phase: Siemens-like 0..4095 (12-bit) stored in uint16
phase_store_max          = 4095;
phase_rescale_slope      = pi/2048;   % radians per stored unit
phase_rescale_intercept  = -pi;       % 0 -> -pi ; 4095 -> ~+pi

% ---------- UIDs ----------
studyUID  = dicomuid;
forUID    = dicomuid;
seriesUID_mag  = arrayfun(@(~) dicomuid, 1:Necho, 'uni', 0);
seriesUID_phs  = arrayfun(@(~) dicomuid, 1:Necho, 'uni', 0);

% ---------- constants ----------
todayStr  = datestr(now,'yyyymmdd');
baseTime  = now;
px = voxelSize_mm(1); py = voxelSize_mm(2); pz = voxelSize_mm(3);
FOVx = Nx*px; FOVy = Ny*py; FOVz = Nz*pz;

% Geometry: IOP axial + IPP at pixel centers
iop = [1 0 0 0 1 0];
x0c = -FOVx/2 + px/2;       % first pixel center (row=1, col=1)
y0c = -FOVy/2 + py/2;

% Common descriptors
manuf   = 'SIEMENS';
model   = 'SyntheticRecon';
softver = 'MATLAB-meFLASH';

% =======================
% MAGNITUDE (one series per echo)
% =======================
fprintf('Exporting MAG to %s\n', outMagDir);
for e = 1:Necho
    thisSeriesUID = seriesUID_mag{e};
    thisTE = TE_ms(e);
    acqTimeStr = datestr(baseTime + (e-1)*seconds(0.1), 'HHMMSS.FFF');

    for z = 1:Nz
        frame = double(img_mag(:, :, z, e));          % [Nx,Ny]
        pix   = uint16(round(frame * mag_scale));     % uint16

        info = struct();
        % SOP / Storage
        info.SOPClassUID                 = '1.2.840.10008.5.1.4.1.1.4';
        info.SOPInstanceUID              = dicomuid;
        info.MediaStorageSOPClassUID     = '1.2.840.10008.5.1.4.1.1.4';
        info.MediaStorageSOPInstanceUID  = info.SOPInstanceUID;

        % IDs
        info.Manufacturer              = manuf;
        info.ManufacturerModelName     = model;
        info.DeviceSerialNumber        = '0000';
        info.SoftwareVersions          = softver;
        info.StudyInstanceUID          = studyUID;
        info.SeriesInstanceUID         = thisSeriesUID;
        info.FrameOfReferenceUID       = forUID;

        % Descriptors
        info.Modality                  = 'MR';
        info.MRAcquisitionType         = '3D';
        info.ImageType                 = 'ORIGINAL\PRIMARY\M\ND';
        info.SeriesDescription         = sprintf('meFLASH_MAG_E%d_TE%.1fms', e, thisTE);
        info.SeriesNumber              = seriesNumberBase + (e-1)*2;   % 800, 802, ...
        info.InstanceNumber            = z;

        % Patient
        info.PatientName               = patientName;
        info.PatientID                 = patientID;
        info.PatientPosition           = 'HFS';
        info.PatientBirthDate          = '';
        info.PatientSex                = '';

        % Times/Dates
        info.StudyDate                 = todayStr;
        info.StudyTime                 = datestr(baseTime,'HHMMSS');
        info.SeriesDate                = todayStr;
        info.SeriesTime                = acqTimeStr;
        info.AcquisitionDate           = todayStr;
        info.AcquisitionTime           = acqTimeStr;
        info.ContentDate               = todayStr;
        info.ContentTime               = acqTimeStr;

        % Geometry
        info.Rows                      = Ny;
        info.Columns                   = Nx;
        info.PixelSpacing              = [py px];        % [row, col]
        info.SliceThickness            = pz;
        info.SpacingBetweenSlices      = pz;
        info.ImageOrientationPatient   = iop;
        zc = -FOVz/2 + (z-0.5)*pz;
        info.ImagePositionPatient      = [x0c; y0c; zc];
        info.SliceLocation             = zc;

        % Acquisition params
        info.RepetitionTime            = TR_ms;
        info.EchoTime                  = thisTE;
        info.EchoNumber                = e;
        info.InPlanePhaseEncodingDirection = 'COL';  % set to 'ROW' if appropriate
        info.FlipAngle                 = 15;         % typical FLASH
        info.NumberOfAverages          = 1;
        info.MagneticFieldStrength     = B0_T;
        info.ImagingFrequency          = 42.576 * B0_T;
        info.PixelBandwidth            = 500;

        % Pixel format
        info.BitsAllocated             = 16;
        info.BitsStored                = 16;
        info.HighBit                   = 15;
        info.PixelRepresentation       = 0;  % unsigned
        info.SamplesPerPixel           = 1;
        info.PhotometricInterpretation = 'MONOCHROME2';
        info.LossyImageCompression     = '00';

        % Real-world scaling
        info.RescaleIntercept          = mag_rescale_intercept;
        info.RescaleSlope              = mag_rescale_slope;
        info.RescaleType               = 'US';
        info.WindowCenter              = 0.25*double(mag_max);
        info.WindowWidth               = 0.5*double(mag_max);

        fname = fullfile(outMagDir, sprintf('IM_MAG_E%02d_SLC%03d.IMA', e, z));
        dicomwrite(pix.', fname, info, 'CreateMode','Copy');  % transpose for DICOM row/col
    end
end

% =======================
% PHASE (one series per echo), Siemens-like storage
% =======================
fprintf('Exporting PHASE to %s\n', outPhsDir);
for e = 1:Necho
    thisSeriesUID = seriesUID_phs{e};
    thisTE = TE_ms(e);
    acqTimeStr = datestr(baseTime + (e-1)*seconds(0.1), 'HHMMSS.FFF');

    for z = 1:Nz
        phi = double(img_phase(:, :, z, e));                     % radians in [-pi, pi)
        stored = round( (phi - phase_rescale_intercept) / phase_rescale_slope ); % (phi + pi)/(pi/2048)
        stored = max(0, min(phase_store_max, stored));           % clamp to [0,4095]
        pix = uint16(stored);                                    % UNSIGNED

        info = struct();
        % SOP / Storage
        info.SOPClassUID                 = '1.2.840.10008.5.1.4.1.1.4';
        info.SOPInstanceUID              = dicomuid;
        info.MediaStorageSOPClassUID     = '1.2.840.10008.5.1.4.1.1.4';
        info.MediaStorageSOPInstanceUID  = info.SOPInstanceUID;

        % IDs
        info.Manufacturer              = manuf;
        info.ManufacturerModelName     = model;
        info.DeviceSerialNumber        = '0000';
        info.SoftwareVersions          = softver;
        info.StudyInstanceUID          = studyUID;
        info.SeriesInstanceUID         = thisSeriesUID;
        info.FrameOfReferenceUID       = forUID;

        % Descriptors
        info.Modality                  = 'MR';
        info.MRAcquisitionType         = '3D';
        info.ImageType                 = 'ORIGINAL\PRIMARY\PHASE\ND';
        info.SeriesDescription         = sprintf('meFLASH_PHASE_E%d_TE%.1fms', e, thisTE);
        info.SeriesNumber              = seriesNumberBase + (e-1)*2 + 1;  % 801, 803, ...
        info.InstanceNumber            = z;

        % Patient
        info.PatientName               = patientName;
        info.PatientID                 = patientID;
        info.PatientPosition           = 'HFS';
        info.PatientBirthDate          = '';
        info.PatientSex                = '';

        % Times/Dates
        info.StudyDate                 = todayStr;
        info.StudyTime                 = datestr(baseTime,'HHMMSS');
        info.SeriesDate                = todayStr;
        info.SeriesTime                = acqTimeStr;
        info.AcquisitionDate           = todayStr;
        info.AcquisitionTime           = acqTimeStr;
        info.ContentDate               = todayStr;
        info.ContentTime               = acqTimeStr;

        % Geometry
        info.Rows                      = Ny;
        info.Columns                   = Nx;
        info.PixelSpacing              = [py px];      % [row, col]
        info.SliceThickness            = pz;
        info.SpacingBetweenSlices      = pz;
        info.ImageOrientationPatient   = iop;
        zc = -FOVz/2 + (z-0.5)*pz;
        info.ImagePositionPatient      = [x0c; y0c; zc];
        info.SliceLocation             = zc;

        % Acquisition params
        info.RepetitionTime            = TR_ms;
        info.EchoTime                  = thisTE;
        info.EchoNumber                = e;
        info.InPlanePhaseEncodingDirection = 'COL';   % set appropriately
        info.FlipAngle                 = 15;
        info.NumberOfAverages          = 1;
        info.MagneticFieldStrength     = B0_T;
        info.ImagingFrequency          = 42.576 * B0_T;
        info.PixelBandwidth            = 500;

        % Pixel format (simulate 12-bit in uint16 container)
        info.BitsAllocated             = 16;
        info.BitsStored                = 12;
        info.HighBit                   = 11;
        info.PixelRepresentation       = 0;    % UNSIGNED for 0..4095
        info.SamplesPerPixel           = 1;
        info.PhotometricInterpretation = 'MONOCHROME2';
        info.LossyImageCompression     = '00';

        % Real-world scaling back to radians
        info.RescaleIntercept          = phase_rescale_intercept; % -pi
        info.RescaleSlope              = phase_rescale_slope;     % pi/2048
        info.RescaleType               = 'PHASE';

        % Viewer window (stored units)
        info.WindowCenter              = 2048;
        info.WindowWidth               = 4096;

        fname = fullfile(outPhsDir, sprintf('IM_PHASE_E%02d_SLC%03d.IMA', e, z));
        dicomwrite(pix.', fname, info, 'CreateMode','Copy');
    end
end

fprintf('DICOM export complete (meFLASH):\n');
fprintf('  Magnitude: %s (%d series, %d slices each)\n', outMagDir, Necho, Nz);
fprintf('  Phase: %s (%d series, %d slices each)\n', outPhsDir, Necho, Nz);
fprintf('  Echo times: ['); fprintf(' %.1f', TE_ms); fprintf(' ] ms\n');
fprintf('  Total files: %d\n', 2 * Necho * Nz);
end
