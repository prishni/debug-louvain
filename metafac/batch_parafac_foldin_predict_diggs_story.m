function batch_predict_diggs_story()
% BATCH_PARAFAC_FOLDIN_PREDICT_DIGG_STORY - runs baseline PARAFAC for digg / comment prediction
%
% model parameters:
%   R_type - which relations are used in the decomposition; modify this
%   setting to run digg prediction ('FullD') or comment prediction ('FullC')
% other configuration:
%   other model / data parameters not specified in this file can be
%   configured in get_digg_option.m

if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);

R_types = {'R2','R23','R123','R234','R1234','R2345','R12345','R23456','R123456', ...
    'R13','R14','R15','R134','R135','R145','R1345', ...
 'FullP','FullD','FullC','FullPCD', ...
 'R4789','R125'};
% R=18,19 - baseline PARAFAC for digg / comment prediction

filename=sprintf('%sblog_%s_v%d-%d_%s%d%s%d',outfilepath,'parafac_foldin_predict_diggs_story',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename '.txt'];
savelog(logfile,logfile);
opts.foldin_sz = 1;

for R=19 %[18 19]
    opts.R_type=R_types{R};
    for K = 4 %[4:4:20]
        opts.K=K;
         opts.pa = 0;
            s = sprintf('predict in R_type=%s K=%d ',opts.R_type,opts.K);
            savelog(logfile,s);
            run_diggs_parafac(opts); % run decomposition by PARAFAC
            predict_diggs_story_parafac(opts); % run fold-in and prediction
    end
end