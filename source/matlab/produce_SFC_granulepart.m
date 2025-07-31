function [successful_completion, error_ca] = produce_SFC_granulepart(varargin)
% Produce part of a 2B-SFC granule, based on the input parameters
%  specified by a 'pd' struct returned from configure_toplevel_IO.
%
% Returns a successful completion boolean flag and a cell array error code.
  
successful_completion = false;   

error_ca = {'#NONE#', 0};  % Default

pd = configure_toplevel_IO;  % Get/set shared filepaths, dirs, control params

% Read in product file data specifications (done here to have access to
%  fill_value parameters for the calculations below):
JSON_fspecs_SFC_fpath = fullfile(pd.ancillary_data_dir, ...
                                 'Sfc_product_filespecs.json');
[SFC_schema, SFC_jdat] = init_nc_schema_from_JSON(JSON_fspecs_SFC_fpath, false);

% Define/initialize output data structure array:
odat = struct;
odat.global_atts = struct;
odat.Geometry = struct;
odat.Geometry.group_atts = struct;
odat.Geometry.group_dims = struct;
odat.Sfc = struct;
%odat.Sfc.group_atts = struct;
odat.Sfc.group_dims = struct;

valid_emis0 = [10 12:16 20:27];   % potential PREFIRE channels for surface emissivity 
% Open input L1B granule data file, read in some useful parameters:
gg_att_C = netcdf.getConstant('NC_GLOBAL');
L1B.ncid = netcdf.open(pd.L1B_rad_fpath, 'NC_NOWRITE');
[L1B.n_dims, ~, L1B.n_globalatts, ~] = netcdf.inq(L1B.ncid);

L1B.G_gid = netcdf.inqNcid(L1B.ncid, 'Geometry');
L1B.R_gid = netcdf.inqNcid(L1B.ncid, 'Radiance');
dimid = netcdf.inqDimID(L1B.ncid, 'atrack');
[~, dims.max_nframes] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'xtrack'); 
[~, dims.nxtrack] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'UTC_parts');
[~, dims.nUTCparts] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'FOV_vertices');
[~, dims.nvertices] = netcdf.inqDim(L1B.ncid, dimid);
dimid = netcdf.inqDimID(L1B.ncid, 'spectral');
[~, dims.nspectral] = netcdf.inqDim(L1B.ncid, dimid);
%sp_subset_ib = 2;  % Spectral subset to be input starts at this (1-based) index 
%dims.nspectral_subset = dims.nspectral-(sp_subset_ib-1);

% Copy some L1B global and group attribute values to the output structure array:
g_atts_to_copy = {'granule_ID', 'spacecraft_ID', 'sensor_ID', ...
                  'ctime_coverage_start_s', 'ctime_coverage_end_s', ...
                  'UTC_coverage_start', 'UTC_coverage_end', ...
                  'orbit_sim_version'};
for ia=1:length(g_atts_to_copy)
   g_att_name = g_atts_to_copy{ia};
   g_att_value = netcdf.getAtt(L1B.ncid, gg_att_C, g_att_name);
   odat.global_atts.(g_att_name) = g_att_value;
end

gp_atts_to_copy = {'image_integration_duration_ms', 'solar_beta_angle_deg', ...
                   'TAI_minus_ctime_at_epoch_s', ...
                   'start_granule_edge_ctime_s', 'end_granule_edge_ctime_s'};
for ia=1:length(gp_atts_to_copy)
   gp_att_name = gp_atts_to_copy{ia};
   odat.Geometry.group_atts.(gp_att_name) = ...
                               netcdf.getAtt(L1B.G_gid, gg_att_C, gp_att_name);
end

% Determine instrument/sensor ID:
sensor_num = str2num(odat.global_atts.sensor_ID(5:6));

% (pd.idxbeg_atrack:pd.idxend_atrack) is the inclusive range/subset to process.
% NOTE that when pd.idxend_atrack == -1, the actual atrack dimension length
%  (determined from the L1B file, since it varies across files) should be used
%  for the end index of the subset.
%
% This means that to process the first N frames, we expect idxbeg_atrack to be 1
% and idxend_atrack to be N.
frame_beg = pd.idxbeg_atrack;
frame_end = pd.idxend_atrack;
if (frame_end == -1)
   % if idxend is -1, that means do the whole file.
   % the inclusive end frame is then the max N frame in the L1B file.
   frame_end = dims.max_nframes;
end
dims.nframes = frame_end-frame_beg+1;

% Read some Geometry group data from the L1B file:
fields_to_read = {'obs_ID', 'latitude', 'longitude'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.G_gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some Geometry group data from the L1B file:
fields_to_read = {'time_UTC_values'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.G_gid, fields_to_read, ...
                                 [1,frame_beg], [dims.nUTCparts,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some Radiance group data from the L1B file:

fields_to_read = {'wavelength', 'idealized_wavelength'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.R_gid, fields_to_read, ...
                  [1,1], [dims.nspectral,dims.nxtrack]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_read = {'spectral_radiance', 'radiance_quality_flag'};
[error_ca, L1B] = read_PREFIRE_granule(L1B, L1B.R_gid, fields_to_read, ...
                  [1,1,frame_beg], [dims.nspectral,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

% Read some meteorology data from the AUX-MET file:

AUXM.ncid = netcdf.open(pd.AUX_MET_fpath, 'NC_NOWRITE');
[AUXM.n_dims, ~, AUXM.n_globalatts, ~] = netcdf.inq(AUXM.ncid);
AUXM.gid = netcdf.inqNcid(AUXM.ncid, 'Aux-Met');
dimid = netcdf.inqDimID(AUXM.ncid, 'zlevels');
[~, dims.nzlevels] = netcdf.inqDim(AUXM.ncid, dimid);
dimid = netcdf.inqDimID(AUXM.ncid, 'n_igbp_classes');
[~, dims.nIGBPc] = netcdf.inqDim(AUXM.ncid, dimid);

fields_to_read = {'surface_pressure', 'skin_temp'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

%fields_to_read = {'VIIRS_surface_type'};
%[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
%                     [1,1,frame_beg], [dims.nIGBPc,dims.nxtrack,dims.nframes]);
%if (error_ca{2} ~= 0)
%   fprintf(2, '%s\n', error_ca{1});
%   return
%end

fields_to_read = {'temp_profile', 'wv_profile', 'o3_profile'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                   [1,1,frame_beg], [dims.nzlevels,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end
g_to_kg = 1.e-3;  % [g -> kg] conversion factor
ppm_to_molrat = 1.e-6;  % [ppm -> molar ratio] conversion factor
molar_mass_O3 = 47.998;  % [g/mol]
molar_mass_air = 28.96*g_to_kg;  % [kg/mol] TO-DO: probably should not be a different value here versus in AUX-MET!
ozone_profile = AUXM.o3_profile*ppm_to_molrat* ...
                     molar_mass_O3/molar_mass_air;  % [ppmv -> gO3/kg_air]

fields_to_read = {'pressure_profile'};
[error_ca, AUXM] = read_PREFIRE_granule(AUXM, AUXM.gid, fields_to_read, ...
                                        [1], [dims.nzlevels]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

%[~, surface_type_xat] = max(AUXM.VIIRS_surface_type, [], 1);

netcdf.close(AUXM.ncid)  % Done with input AUX-MET file

% Read some cloud mask data from the 2B-MSK file:

MSK.ncid = netcdf.open(pd.L2B_msk_fpath, 'NC_NOWRITE');
[MSK.n_dims, ~, MSK.n_globalatts, ~] = netcdf.inq(MSK.ncid);
MSK.gid = netcdf.inqNcid(MSK.ncid, 'Msk');

fields_to_read = {'cldmask_probability'};
[error_ca, MSK] = read_PREFIRE_granule(MSK, MSK.gid, fields_to_read, ...
                                   [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(MSK.ncid)  % Done with input 2B-MSK file

if sensor_num == 1
   data_pri_emis = load(fullfile(pd.ancillary_data_dir, ...
                                 'emis0deg_allyear_60-90N_0-360E_TIRS1.mat'));
elseif sensor_num ==2
   data_pri_emis = load(fullfile(pd.ancillary_data_dir, ...
                                 'emis0deg_allyear_60-90N_0-360E_TIRS2.mat'));
end
prior_emis_data = data_pri_emis.emis0deg_prefire;

% Read useful data from SRF file:
SRF_fpath = fullfile(pd.ancillary_data_dir, ...
            sprintf('PREFIRE_TIRS%1d_%s.nc', sensor_num, pd.SRF_disambig_str));
ncid = netcdf.open(SRF_fpath, 'NC_NOWRITE');
fields_to_read = {'channel_wavenum1', 'channel_wavenum2', ...
                  'channel_wavelen1', 'channel_wavelen2', 'NEDR'};
for iv=1:length(fields_to_read)
   v_name = fields_to_read{iv};
   varid = netcdf.inqVarID(ncid, v_name);
   SRF.([v_name, '_T']) = netcdf.getVar(ncid, varid)';  % C-indexing
end
fields_to_read = {'srf', 'wavelen', 'detector_bitflags'};
for iv=1:length(fields_to_read)
   v_name = fields_to_read{iv};
   varid = netcdf.inqVarID(ncid, v_name);
   SRF.(v_name) = netcdf.getVar(ncid, varid);  % C-indexing
end
SRF.im_version = netcdf.getAtt(ncid, gg_att_C, 'instrument_model_version');
netcdf.close(ncid)  % Done with input SRF file

% Read useful data from MODTRAN standard profiles file:
MODTRAN_default_profile = load(fullfile(pd.ancillary_data_dir, ...
                                        'Modtran_standard_profiles.mat'));

% Get the "stock" PCRTM wavenumbers and emissivities in the case of sensor id2
%  (number of channels 740, wavenumber range 50.38375:2759.88625)
fid = fopen(fullfile(pd.ancillary_data_dir, 'EMISSIVITY_id2_copy'), 'r');
c = textscan(fid, '%s', 'Delimiter', '\n');
emis_id2.cell_array_of_lines = c{1};
nhd = 18;  % Number of "header" lines
nch = 740;  % Number of channel entries in the file
emis_id2.wn_new = zeros(nch,1);
for i=1:nch
   data = str2num(emis_id2.cell_array_of_lines{i+nhd});
   if i == 1
      nvl = length(data)-1;
      emis_id2.emis = zeros(nvl,nch);
   end
   emis_id2.wn_new(i) = data(1);
   emis_id2.emis(:,i) = data(2:end);
end
fclose(fid);

% Find all channel indices with idealized wavelength < 5 um:
idw = zeros(dims.nspectral, dims.nxtrack);
n_idw = zeros(dims.nxtrack);
for ix=1:dims.nxtrack
   itmp = find(L1B.idealized_wavelength(:,ix) < 5);
   n_idw(ix) = length(itmp);
   idw(1:n_idw(ix),ix) = itmp;
end

%-- Some constraints about which FOVs a retrieval will be attempted for:

% Ch 14 is good on both SAT1 and SAT2; for certain quality-check operations,
%  this channel is selected in order to yield a boolean array of the correct
%  shape that reflects the FOV (i.e., cross-track + along-track) dependence of
%  a field.
qc_ch = 14;

% These are per FOV (dimensioned {dims.nxtrack,dims.nframes}):
% (note that AUX-MET info is created for every FOV present in a 1B-RAD granule,
%  and so does not have one of these boolean masks)
lat_okay = ( abs(L1B.latitude) > 60. );
cmsk_okay = ( MSK.cldmask_probability < 0.4 );
cmsk_better = ( MSK.cldmask_probability < 0.1 );
L1B_FOV_okay = squeeze(L1B.radiance_quality_flag(qc_ch,:,:) <= 1);
attempt_retrieval = ( lat_okay & cmsk_okay & L1B_FOV_okay );

n_emis_gt_maxth = zeros(dims.nxtrack,dims.nframes);
n_emis_gt_unity = zeros(dims.nxtrack,dims.nframes);
n_emis_lt_minth = zeros(dims.nxtrack,dims.nframes);
iter = zeros(dims.nxtrack,dims.nframes);
ec = zeros(dims.nxtrack,dims.nframes);
emis = NaN(dims.nspectral,dims.nxtrack,dims.nframes);
emis_unc = NaN(dims.nspectral,dims.nxtrack,dims.nframes);
for ia=1:dims.nframes
   for ix=1:dims.nxtrack
%%%% each nxtrack is a scene
      footprint_id = ix;
      mask = squeeze(SRF.detector_bitflags(footprint_id,:));
      idx = find(mask(valid_emis0) < 1);  % Only use good, unmasked channels
      clear valid_emis; valid_emis = valid_emis0(idx);  
      clear valid_rad; valid_rad = valid_emis;

      lat = L1B.latitude(ix,ia);
      if attempt_retrieval(ix,ia)

%         fprintf('Starting xtrack=%d, atrack=%d, lat= %f ...', ix, ...
%                 ia+frame_beg-1, lat);

         [wv_emis, wv_rad, wn_emis, wn_rad, x, y, S_aposteriori, x_op, y_op, ...
          S_op, K_op, dgf, converged, iter(ix,ia), chi2, pvalue, c2, ...
          ec(ix,ia)] = OEemis(pd.PCRTM_dir, pd.ancillary_data_dir, ...
                   ix, valid_rad, valid_emis, ...
                   L1B.spectral_radiance(:,ix,ia), lat, ...
                   L1B.longitude(ix,ia), 0, ...
                   AUXM.pressure_profile', AUXM.surface_pressure(ix,ia), ...
                   AUXM.wv_profile(:,ix,ia)', ozone_profile(:,ix,ia)', ...
                   AUXM.temp_profile(:,ix,ia)', AUXM.skin_temp(ix,ia), ...
                   prior_emis_data, SRF, MODTRAN_default_profile, emis_id2, ...
                   ec(ix,ia));

         emis(:,ix,ia) = expandEmis(x_op, SRF, ix, valid_emis);
         emis_unc(:,ix,ia) = expandEmis(sqrt(diag(S_op)), SRF, ix, valid_emis);

         % Set any emissivity values at wavelengths < 5 um to missing values:
         emis(idw(1:n_idw(ix),ix),ix,ia) = NaN;
         emis_unc(idw(1:n_idw(ix),ix),ix,ia) = NaN;

         % Make a record of each FOV that has any of the following "special"
         %  spectral emissivity value ranges:
         n_emis_gt_unity(ix,ia) = sum(emis(:,ix,ia) > 1.);
         n_emis_gt_maxth(ix,ia) = sum(emis(:,ix,ia) > 1.1);
         n_emis_lt_minth(ix,ia) = sum(emis(:,ix,ia) < 0.5);
      end
   
   end
end

%=== NOTE: The order of the following quality checks/assignments is important:

sfc_quality_flag = zeros(dims.nxtrack, dims.nframes, 'int8');
qc_bitflags = zeros(dims.nxtrack, dims.nframes, 'uint16');

% Set any emissivity channels that have NaN values to the proper output
%  FillValue.  If all channels are NaN for a given FOV, set the quality flag
%  to its appropriate FillValue as well.
NaN_bool = isnan(emis);
emis(NaN_bool) = SFC_jdat.Sfc.sfc_spectral_emis.fill_value;
emis_unc(NaN_bool) = SFC_jdat.Sfc.sfc_spectral_emis_unc.fill_value;
sfc_quality_flag(squeeze(all(NaN_bool, 1))) = ...
                                      SFC_jdat.Sfc.sfc_quality_flag.fill_value;

% Set some "why retrieval was not attempted" bitflags:
idx = find(~lat_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 0+1, 'uint16');
idx = find(~L1B_FOV_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 1+1, 'uint16');
idx = find(~cmsk_okay);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 2+1, 'uint16');

% Set bitflag indicating where retrievals were performed with
%  cldmsk_probability < 0.1:
idx = find(cmsk_better);
qc_bitflags(idx) = bitset(qc_bitflags(idx), 10+1, 'uint16');

% Fill OE_iterations output field; where iterations=0, set to FillValue:
odat.Sfc.OE_iterations = int8(iter);
odat.Sfc.OE_iterations(odat.Sfc.OE_iterations == 0) = ...
                                    int8(SFC_jdat.Sfc.OE_iterations.fill_value);

% Record some iteration error/warning codes in bitflags field:
idx = find(ec == 2);  % negative convergence criterion at last iteration
qc_bitflags(idx) = bitset(qc_bitflags(idx), 3+1, 'uint16');
idx = find(ec == 3);  % zero degrees of freedom at last iteration
qc_bitflags(idx) = bitset(qc_bitflags(idx), 4+1, 'uint16');

%-- Set FOVs with 1+ channel(s) that exceed set threshold emissivity values
%    to FillValue, setting flags and bitflags accordingly:

[ee_row, ee_col] = find(n_emis_gt_maxth > 0);
emis(:, ee_row, ee_col) = SFC_jdat.Sfc.sfc_spectral_emis.fill_value;
emis_unc(:, ee_row, ee_col) = SFC_jdat.Sfc.sfc_spectral_emis_unc.fill_value;
sfc_quality_flag(ee_row, ee_col) = SFC_jdat.Sfc.sfc_quality_flag.fill_value;
b_idx = (n_emis_gt_maxth > 0) & (n_emis_gt_maxth <= 2);
qc_bitflags(b_idx) = bitset(qc_bitflags(b_idx), 5+1, 'uint16');
b_idx = n_emis_gt_maxth > 2;
qc_bitflags(b_idx) = bitset(qc_bitflags(b_idx), 6+1, 'uint16');

[ee_row, ee_col] = find(n_emis_lt_minth > 0);
emis(:, ee_row, ee_col) = SFC_jdat.Sfc.sfc_spectral_emis.fill_value;
emis_unc(:, ee_row, ee_col) = SFC_jdat.Sfc.sfc_spectral_emis_unc.fill_value;
sfc_quality_flag(ee_row, ee_col) = SFC_jdat.Sfc.sfc_quality_flag.fill_value;
b_idx = (n_emis_lt_minth > 0) & (n_emis_lt_minth <= 2);
qc_bitflags(b_idx) = bitset(qc_bitflags(b_idx), 7+1, 'uint16');
b_idx = n_emis_lt_minth > 2;
qc_bitflags(b_idx) = bitset(qc_bitflags(b_idx), 8+1, 'uint16');

% Set flag and bitflag for FOVs with 1+ channel(s) that exceed emissivity = 1:
b_idx = n_emis_gt_unity > 0;
qc_bitflags(b_idx) = bitset(qc_bitflags(b_idx), 9+1, 'uint16');
sfc_quality_flag(b_idx) = 1;

% Fill sfc_quality_flag and sfc_qc_bitflags output fields;
odat.Sfc.sfc_quality_flag = int8(sfc_quality_flag);
odat.Sfc.sfc_qc_bitflags = qc_bitflags;

% Load any additional data (including dimension lengths) into output structure
%  array(s):

odat.global_atts.full_versionID = pd.product_fullversion;

line_parts = strsplit(pd.product_fullversion, '_');
odat.global_atts.archival_versionID = strrep(line_parts{1}, 'R', '');

tmp_fpath = fullfile(pd.top_path, 'dist', ...
                     sprintf('prdgit_version_m%d.txt', pd.proc_mode));
fid = fopen(tmp_fpath, 'rt');
tmp_line = fgetl(fid);
fclose(fid);
idx = strfind(tmp_line, '(');
odat.global_atts.provenance = sprintf('%s%s (%s', ...
                   extractBefore(tmp_line, idx(1)), pd.product_fullversion, ...
                   extractAfter(tmp_line, idx(1)));

iftmp = {pd.L1B_rad_fpath, pd.AUX_MET_fpath, pd.L2B_msk_fpath};
if (pd.proc_mode == 2)
   iftmp{4} = pd.L2B_atm_fpath;
end
for ia=1:length(iftmp)
   [~, name, ext] = fileparts(iftmp{ia});
   intmp{ia} = [name ext];
end
odat.global_atts.input_product_files = strjoin(intmp, ', ');

odat.global_atts.SRF_NEdR_version = SRF.im_version;

tmp_fpath = fullfile(pd.top_path, 'VERSION.txt');
fid = fopen(tmp_fpath, 'rt');
odat.global_atts.processing_algorithmID = fgetl(fid);
fclose(fid);

odat.Geometry.group_dims.atrack = dims.nframes;
odat.Geometry.group_dims.xtrack = dims.nxtrack;
odat.Geometry.group_dims.UTC_parts = dims.nUTCparts;
odat.Geometry.group_dims.FOV_vertices = dims.nvertices;

fields_to_copy = {'obs_ID', 'latitude', 'longitude', 'time_UTC_values'};
for i=1:length(fields_to_copy)
    odat.Geometry.(fields_to_copy{i}) = L1B.(fields_to_copy{i});
end

fields_to_readcopy = {'ctime', 'ctime_minus_UTC', 'subsat_latitude', ...
                      'subsat_longitude', 'sat_altitude', ...
                      'sat_solar_illumination_flag', 'orbit_phase_metric', ...
                      'satellite_pass_type'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                                 [frame_beg], [dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_readcopy = {'vertex_latitude', 'vertex_longitude', ...
                      'maxintgz_verts_lat', 'maxintgz_verts_lon'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                  [1,1,frame_beg], [dims.nvertices,dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

fields_to_readcopy = {'land_fraction', ...
                   'elevation', 'elevation_stdev', 'viewing_zenith_angle', ...
                   'viewing_azimuth_angle', 'solar_zenith_angle', ...
                   'solar_azimuth_angle', 'solar_distance', ...
                   'geoloc_quality_bitflags'};
[error_ca, odat.Geometry] = read_PREFIRE_granule(odat.Geometry, L1B.G_gid, ...
                                 fields_to_readcopy, ...
                                 [1,frame_beg], [dims.nxtrack,dims.nframes]);
if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

netcdf.close(L1B.ncid);  % Done with input L1B file

odat.Sfc.group_dims.atrack = dims.nframes;
odat.Sfc.group_dims.xtrack = dims.nxtrack;
odat.Sfc.group_dims.spectral = dims.nspectral;

odat.Sfc.idealized_wavelength = L1B.idealized_wavelength;
odat.Sfc.wavelength = L1B.wavelength;

odat.Sfc.sfc_spectral_emis = emis;
odat.Sfc.sfc_spectral_emis_unc = emis_unc;

now_UTCn_DT = datetime('now', 'TimeZone', 'UTC', 'Format', ...
                      'yyyy-MM-dd''T''HH:mm:ss.SSSSSS');  % faux-UTC (no leap-s)
odat.global_atts.UTC_of_file_creation = string(now_UTCn_DT);

% UTC shape is (7, n), where the 7th row is millisec.
% change to a UTC array shaped (6, n) where the 6th row is
% fractional seconds.
UTC_granule_start = L1B.time_UTC_values(1:6,1);  % int16 array

% Construct output filepath, inserting a suffix to describe the
% 'granule part' (as an inclusive atrack range.)
fn_fmtstr = 'raw-PREFIRE_SAT%1d_2B-SFC_%s_%04d%02d%02d%02d%02d%02d_%s';
output_fn = sprintf(fn_fmtstr, sensor_num, pd.product_fullversion, ...
                    UTC_granule_start(1:6,1), odat.global_atts.granule_ID);
suffix_fmtstr = '-%s_%05d_%05d_of_%05df.nc';

% Output frame range (in filename) should be the original inclusive
%  0-index values.
if pd.idxend_atrack == -1
   iend = frame_end-1;  % 0-based index
else
   iend = pd.idxend_atrack-1;  % 0-based index
end
part_suffix = sprintf(suffix_fmtstr, 'atrack', pd.idxbeg_atrack-1, iend, ...
                      dims.max_nframes);
output_fpath = fullfile(pd.output_dir, [output_fn part_suffix]);

[~, name, ~] = fileparts(output_fpath);
odat.global_atts.file_name = name;

% Prepare for NetCDF-format file write, then define and write the file:
[m_ncs, vars_to_write] = modify_nc_schema(SFC_schema, odat);
if isfile(output_fpath)
   delete(output_fpath);
end
ncwriteschema(output_fpath, m_ncs);
write_predef_nc_vars(output_fpath, vars_to_write);

if (error_ca{2} ~= 0)
   fprintf(2, '%s\n', error_ca{1});
   return
end

successful_completion = true;

end
