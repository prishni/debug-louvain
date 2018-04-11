function run_diggs_parafac(opts)

% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s] = get_digg_option(opts);

filename=sprintf('%slog_%s_v%d-%d_%s%d%s%d',outfilepath,'run_diggs_parafact',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename 'K' num2str(K) R_type 'pa' num2str(pa) '.txt'];
savelog(logfile,logfile);
file_prefix=sprintf('%sresPARAFAC_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);

tmin=TMIN; tmax=TMAX;
for ti = tmin:tmax
    if ti<tmin,continue;end  
    % load data tensors
%     filename=sprintf('%ss%s%s%d_v%d.mat',infilepath,R_type,stream_s,ti,data_ver);
    switch R_type
        case {'FullD'}
        filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        X = sFullD; clear sFullD;
        case {'FullC'}
        filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullC=sptensor(sFullC);% get sD (user,story) via digg
        X = sFullC; clear sFullC;
        case {'FullDCP'}
        filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        X = sFullD; clear sFullD;
        filename=sprintf('%ssFullP%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullP=sptensor(sFullP);% get sD (user,story) via digg
        X = X+sFullP; clear sFullP;
        filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullC=sptensor(sFullC);% get sD (user,story) via digg
        X = X+sFullC; clear sFullC;
    end
    tic;
    P=parafac_als(X,K);
    S = P.lambda;
    U = P.U;
    cvtime = toc;

    decom_res_file = [file_prefix 'K' num2str(K) R_type 't' num2str(ti) '.mat'];    
    save(decom_res_file,'S','U','cvtime');  
    fprintf('\nsave %s',decom_res_file);
    savelog(logfile, sprintf('save %s',decom_res_file));        
end