function ssfp2dcm(prefix, img_mag, img_phase, TE0_ms, TR_ms, voxelSize_mm, B0_T, patientName, patientID, seriesNumberBase)
% Export SSFP magnitude and phase images to DICOM with corrected orientation
% and Siemens-like phase scaling.
%
% prefix         : output folder prefix (created if missing), e.g. 'out/tr8_n5'
% img_mag/phase  : single, size [Nx,Ny,Nz,Necho]; phase in radians (-pi..pi)
% TE0_ms/TR_ms   : ms (TE at center of readout)
% voxelSize_mm   : [dx dy dz]
% B0_T           : Tesla, e.g. 3
% patientName/ID : strings for DICOM tags (can be placeholders)
% seriesNumberBase : integer (e.g., 800 for mag series, phase will use base+1)
%
% Key choices:
% - Phase stored UNSIGNED 0..4095 (uint16) with RescaleSlope=pi/2048, RescaleIntercept=-pi
%   so viewers map back to [-pi, +pi).
% - ImageType set to Siemens-like strings ('...\\M\\ND' and '...\\PHASE\\ND').
% - EchoNumber = e for each F-state series.
% - Geometry: IOP = [1 0 0 0 1 0], IPP at pixel centers.

arguments
    prefix (1,:) char
    img_mag single
    img_phase single
    TE0_ms (1,1) double
    TR_ms (1,1) double
    voxelSize_mm (1,3) double
    B0_T (1,1) double = 3
    patientName (1,:) char = 'Anon^SSFP'
    patientID   (1,:) char = '000000'
    seriesNumberBase (1,1) double = 0
end

% ---------- shapes ----------
sz = size(img_mag); if numel(sz)==3, sz(4)=1; end
[Nx,Ny,Nz,Necho] = deal(sz(1),sz(2),sz(3),sz(4));
assert(isequal(size(img_phase), [Nx,Ny,Nz,Necho]), 'mag/phase shape mismatch');

%% ---------- Apply orientation corrections ONCE before export ----------
% Flip L-R (dim 2) and slice order (dim 3) to match typical Siemens axial
fprintf('Applying orientation corrections...\n');
img_mag_corrected   = zeros(size(img_mag),   'single');
img_phase_corrected = zeros(size(img_phase), 'single');

for e = 1:Necho
    tmpm = flip(flip(img_mag(:,:,:,e), 2), 3);
    tmpp = flip(flip(img_phase(:,:,:,e), 2), 3);
    img_mag_corrected(:,:,:,e)   = tmpm;
    img_phase_corrected(:,:,:,e) = tmpp;
end

img_mag   = img_mag_corrected;
img_phase = img_phase_corrected;

% ---------- output dirs ----------
outMagDir = fullfile(prefix,'DICOM_MAG');   if ~exist(outMagDir,'dir'),   mkdir(outMagDir);   end
outPhsDir = fullfile(prefix,'DICOM_PHASE'); if ~exist(outPhsDir,'dir'), mkdir(outPhsDir); end

% ---------- scaling ----------
% Magnitude stored as uint16 with linear rescale back to original range
mag_max   = max(img_mag(:));  if mag_max<=0, mag_max = 1; end
mag_scale = 65535 / double(mag_max);    % store as uint16
mag_rescale_slope   = 1 / mag_scale;
mag_rescale_intercept = 0;

% Phase stored Siemens-like: 0..4095 in uint16, maps back via slope/intercept
phase_store_max      = 4095;
phase_rescale_slope  = pi/2048;   % radians per stored unit
phase_rescale_intercept = -pi;    % so 0 -> -pi, 4095 -> ~+pi

% ---------- UIDs ----------
studyUID  = dicomuid;
forUID    = dicomuid;
seriesUID_mag  = arrayfun(@(~) dicomuid, 1:Necho, 'uni', 0);
seriesUID_phs  = arrayfun(@(~) dicomuid, 1:Necho, 'uni', 0);

% ---------- constants ----------
todayStr  = datestr(now,'yyyymmdd');
baseTime  = now;
px = voxelSize_mm(1); py = voxelSize_mm(2); pz = voxelSize_mm(3);

% DICOM geometry
% IOP for standard axial: rows L->R (x), cols P->A (y), normal S (z)
iop = [1 0 0 0 1 0];

% Helper for IPP at pixel centers for slice z (1-based)
% First pixel center (row=1,col=1) of each slice:
% x0 = -FOVx/2 + px/2 ; y0 = -FOVy/2 + py/2 ; z0 = -FOVz/2 + (z-0.5)*pz
FOVx = Nx*px; FOVy = Ny*py; FOVz = Nz*pz;
x0c  = -FOVx/2 + px/2;
y0c  = -FOVy/2 + py/2;

% Common tags
manuf   = 'SIEMENS';
model   = 'SyntheticRecon';
softver = 'MATLAB-SSFP';

% =======================
% Write MAGNITUDE (one series per echo)
% =======================
fprintf('Exporting magnitude DICOM images to %s\n', outMagDir);
for e = 1:Necho
    fprintf('  Magnitude: Echo %d/%d\n', e, Necho);
    thisSeriesUID = seriesUID_mag{e};
    thisTE = TE0_ms;
    acqTimeStr = datestr(baseTime + (e-1)*seconds(0.1), 'HHMMSS.FFF');

    for z = 1:Nz
        frame = double(img_mag(:, :, z, e));         % [Nx,Ny] MATLAB
        pix   = uint16(round(frame * mag_scale));    % uint16

        info = struct();
        % SOP / Storage
        info.SOPClassUID                 = '1.2.840.10008.5.1.4.1.1.4';
        info.SOPInstanceUID              = dicomuid;
        info.MediaStorageSOPClassUID     = '1.2.840.10008.5.1.4.1.1.4';
        info.MediaStorageSOPInstanceUID  = info.SOPInstanceUID;

        % Identification
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
        info.SeriesDescription         = sprintf('SSFP_MAG_F%d_TE%.1fms', e-1, thisTE);
        info.SeriesNumber              = seriesNumberBase + (e-1)*2; % 800,802,804,...
        info.InstanceNumber            = z;

        % Patient
        info.PatientName               = patientName;
        info.PatientID                 = patientID;
        info.PatientPosition           = 'HFS';
        info.PatientBirthDate          = '';
        info.PatientSex                = '';

        % Dates/Times
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
        info.PixelSpacing              = [py px];     % [row, col]
        info.SliceThickness            = pz;
        info.SpacingBetweenSlices      = pz;
        info.ImageOrientationPatient   = iop;
        zc = -FOVz/2 + (z-0.5)*pz;
        info.ImagePositionPatient      = [x0c; y0c; zc];
        info.SliceLocation             = zc;

        % Acquisition params
        info.RepetitionTime            = TR_ms;
        info.EchoTime                  = thisTE;
        info.EchoNumber                = e;          % 1..Necho
        info.InPlanePhaseEncodingDirection = 'COL';  % set to 'ROW' if appropriate
        info.FlipAngle                 = 8;
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
        info.WindowCenter              = 0.25*double(mag_max);   % viewer hint
        info.WindowWidth               = 0.5*double(mag_max);    % viewer hint

        fname = fullfile(outMagDir, sprintf('IM_MAG_F%02d_SLC%03d.IMA', e-1, z));
        % Transpose to swap row/col into DICOM order; use .' for real transpose
        dicomwrite(pix.', fname, info, 'CreateMode','Copy');
    end
end

% =======================
% Write PHASE (one series per echo) Siemens-like
% =======================
fprintf('Exporting phase DICOM images to %s\n', outPhsDir);
for e = 1:Necho
    fprintf('  Phase: Echo %d/%d\n', e, Necho);
    thisSeriesUID = seriesUID_phs{e};
    thisTE = TE0_ms;
    acqTimeStr = datestr(baseTime + (e-1)*seconds(0.1), 'HHMMSS.FFF');

    for z = 1:Nz
        phi = double(img_phase(:, :, z, e));       % radians in [-pi, pi)
        stored = round( (phi - phase_rescale_intercept) / phase_rescale_slope ); % (phi + pi)/(pi/2048)
        stored = max(0, min(phase_store_max, stored));   % clamp to [0,4095]
        pix = uint16(stored);                      % UNSIGNED

        info = struct();
        % SOP / Storage
        info.SOPClassUID                 = '1.2.840.10008.5.1.4.1.1.4';
        info.SOPInstanceUID              = dicomuid;
        info.MediaStorageSOPClassUID     = '1.2.840.10008.5.1.4.1.1.4';
        info.MediaStorageSOPInstanceUID  = info.SOPInstanceUID;

        % Identification
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
        info.SeriesDescription         = sprintf('SSFP_PHASE_F%d_TE%.1fms', e-1, thisTE);
        info.SeriesNumber              = seriesNumberBase + (e-1)*2 + 1; % 801,803,805,...
        info.InstanceNumber            = z;

        % Patient
        info.PatientName               = patientName;
        info.PatientID                 = patientID;
        info.PatientPosition           = 'HFS';
        info.PatientBirthDate          = '';
        info.PatientSex                = '';

        % Dates/Times
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
        info.PixelSpacing              = [py px];     % [row, col]
        info.SliceThickness            = pz;
        info.SpacingBetweenSlices      = pz;
        info.ImageOrientationPatient   = iop;
        zc = -FOVz/2 + (z-0.5)*pz;
        info.ImagePositionPatient      = [x0c; y0c; zc];
        info.SliceLocation             = zc;

        % Acquisition params
        info.RepetitionTime            = TR_ms;
        info.EchoTime                  = thisTE;
        info.EchoNumber                = e;          % 1..Necho
        info.InPlanePhaseEncodingDirection = 'COL';  % set appropriately
        info.FlipAngle                 = 8;
        info.NumberOfAverages          = 1;
        info.MagneticFieldStrength     = B0_T;
        info.ImagingFrequency          = 42.576 * B0_T;
        info.PixelBandwidth            = 500;

        % Pixel format (simulate 12-bit within uint16 container)
        info.BitsAllocated             = 16;
        info.BitsStored                = 12;   % Siemens often uses 12-bit for phase
        info.HighBit                   = 11;
        info.PixelRepresentation       = 0;    % UNSIGNED (critical for 0..4095)
        info.SamplesPerPixel           = 1;
        info.PhotometricInterpretation = 'MONOCHROME2';
        info.LossyImageCompression     = '00';

        % Real-world scaling back to radians
        info.RescaleIntercept          = phase_rescale_intercept; % -pi
        info.RescaleSlope              = phase_rescale_slope;     % pi/2048
        info.RescaleType               = 'PHASE';

        % Viewer window (in stored units)
        info.WindowCenter              = 2048;    % mid of 0..4095
        info.WindowWidth               = 4096;

        fname = fullfile(outPhsDir, sprintf('IM_PHASE_F%02d_SLC%03d.IMA', e-1, z));
        dicomwrite(pix.', fname, info, 'CreateMode','Copy');
    end
end

fprintf('DICOM export complete:\n');
fprintf('  Magnitude: %s (%d series, %d slices each)\n', outMagDir, Necho, Nz);
fprintf('  Phase: %s (%d series, %d slices each)\n', outPhsDir, Necho, Nz);
fprintf('  Total files: %d\n', 2 * Necho * Nz);
end
