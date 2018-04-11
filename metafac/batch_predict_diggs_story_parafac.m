function batch_run_diggs_parafac()
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s] = get_digg_option(opts);

% my_startup_20090127 % change filenames

R_types = {'FullP','FullD','FullC','FullDCP'};

filename=sprintf('%sblog_%s_v%d-%d_%s%d%s%d',outfilepath,'predict_diggs_story_parafac',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename '.txt'];
savelog(logfile,logfile);

for R=4
    opts.R_type=R_types{R};
    for K = 4:4:20
        opts.K=K;
        s = sprintf('running R_type=%s K=%d ...',opts.R_type,opts.K);
        savelog(logfile,s);
        predict_diggs_story_parafac(opts);
    end
end