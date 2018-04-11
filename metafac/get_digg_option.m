function [K,pa,TMIN,TMAX,ndt,R_type, ...
    infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s ...
    topN,topP,foldin_sz] = get_digg_option(opts)
% GET_DIGG_OPTION - configure the experiment setting, including model /
% data parameters.
%
%   model parameters: 
%     K - number of factors to be decomposed 
%     pa - weight for prior model
%     R_type - which relations are used in the decomposition
%     TMIN / TMAX - the indices of start and end time of the input data
%     ndt - data interval (not used)
%   data parameters:
%     infilepath0 - input filepath, consisting of files for IDs
%     corresponding to Digg objects (users, tags, etc.)
%     infilepath - input filepath, consisting of tensor data files
%     outfilepath - output filepath
%     stream - the time window for segmenting data; available options {'month','week','3day'}
%     data_ver / res_ver - the version number of data and result
%
% Experiemnt data:
%   The testing data can be downloaded at
%   http://www.public.asu.edu/~ylin56/kdd09sup.html
%   Modify the input filepaths (infilepath0 / infilepath) in the
%   init_option() function below to run experiments on the downloaded
%   dataset.

% Author: Yu-Ru Lin <yu-ru.lin@asu.edu>, Feb 2009

if nargin == 0,opts = struct();end
opts = init_option(opts);
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_option(opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,pa,TMIN,TMAX,ndt,R_type, ...
    infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s ...
    topN,topP,foldin_sz] = get_option(opts)

% model params
K=opts.K; pa = opts.pa; 
TMIN=opts.TMIN; TMAX = opts.TMAX; ndt=opts.ndt;
R_type = opts.R_type;

% data params
infilepath0=opts.infilepath0;
infilepath=opts.infilepath;
outfilepath=opts.outfilepath;
data_ver=opts.data_ver;
res_ver=opts.res_ver;
stream=opts.stream;
if strcmp(stream,'month')==1,stream_s='m';
elseif strcmp(stream,'week')==1,stream_s='w';
elseif strcmp(stream,'3day')==1,stream_s='3d';
end
topN = opts.topN;
topP = opts.topP;
foldin_sz = opts.foldin_sz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function opts = init_option(opts)
if ~isfield(opts, 'K'), opts.K = 12; end
if ~isfield(opts, 'pa'), opts.pa = .2; end
if ~isfield(opts, 'TMIN'), opts.TMIN = 1; end
if ~isfield(opts, 'TMAX'), opts.TMAX = 9; end
if ~isfield(opts, 'ndt'), opts.ndt = 1; end
if ~isfield(opts, 'R_type'), opts.R_type = 'FullDCP'; end
%if ~isfield(opts, 'R_type'), opts.R_type = 'freqC'; end

% modify opts.infilepath0 and opts.infilepath to specify input paths
if ~isfield(opts, 'infilepath0'), opts.infilepath0 = '../data/texts/'; end
if ~isfield(opts, 'infilepath'), opts.infilepath = '../data/sptensors/'; end
% modify opts.outfilepath to specify output path
if ~isfield(opts, 'outfilepath'), opts.outfilepath = 'output_full/'; end
if ~isfield(opts, 'data_ver'), opts.data_ver = 1; end
if ~isfield(opts, 'res_ver'), opts.res_ver = 1; end
if ~isfield(opts, 'stream'), opts.stream = '3day'; end

% disp
if ~isfield(opts, 'topN'), opts.topN = 5; end

% prediction
if ~isfield(opts, 'topP'), opts.topP = 10; end
if ~isfield(opts, 'foldin_sz'), opts.foldin_sz = 1; end
% if foldin_sz, compute U_s = p(s|z); else U_s = p(z|s)

rand('seed',21);
randn('seed',2100);