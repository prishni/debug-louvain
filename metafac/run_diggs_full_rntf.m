function run_diggs_full_rntf(opts)
% run hypergraph factorization on Diggs FULL TENSOR data

% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s] = get_digg_option(opts);
% pa=0; % note: baseline should not have informative prior

filename=sprintf('%slog_%s_v%d-%d_%s%d%s%d',outfilepath,'run_diggs_rntf',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename 'K' num2str(K) R_type 'pa' num2str(pa) '.txt'];
lfile = fopen(logfile, 'w');
fprintf(lfile,logfile);
file_prefix=sprintf('%sres_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);

tstart=TMIN; tmax=TMAX;
for ti = 1:tmax
    if ti<tstart,continue;end  
    % load data tensors
    switch(R_type)
        case {'FullP'} % user,story,word,topic, via submit
            filename=sprintf('%ssFullP%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sptensor(sFullP); clear sFullP; % get sP (user,story) via submit
        case {'FullD'} % via digg
            filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sptensor(sFullD); clear sFullD; % get sD (user,story) via digg
        case {'FullC'} % via comment
            filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sptensor(sFullC); clear sFullC; % get sC (user,story,comment) 
        case {'FullDCP'} % via comment
            filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sptensor(sFullD); clear sFullD; % get sD (user,story) via digg
            filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sT+sptensor(sFullC); clear sFullC; % get sC (user,story,comment) 
            filename=sprintf('%ssFullP%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
            load(filename); sT=sT+sptensor(sFullP); clear sFullP; % get sP (user,story) via submit
    end

    all_dims = size(sT);
%     all_dims = [nU,nS,nW,nJ];
    fprintf('\nt=%d',ti);
    if ti>tstart,prior = {S,U,pa};
    else prior={}; end
    
    decom_res_file = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 't' num2str(ti) '.mat'];
    if exist(decom_res_file,'file')>0
        load(decom_res_file,'S','U');  
        fprintf('\n%s exist',decom_res_file);
        fprintf(lfile, sprintf('%s exist',decom_res_file));
    else
        Vdims = all_dims;
        RG = {{sT,1:4,1}};
        [S,U,iters,cvtime,ll]=rntf(RG,Vdims,K,prior);
        save(decom_res_file,'S','U','iters','cvtime','ll');  
        fprintf('\nsave %s',decom_res_file);
        fprintf(lfile, sprintf('save %s',decom_res_file));    
    end    
end