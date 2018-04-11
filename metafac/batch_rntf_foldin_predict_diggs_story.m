function batch_predict_diggs_story()
% BATCH_RNTF_FOLDIN_PREDICT_DIGG_STORY - runs MetaFac and baseline MWA for digg / comment prediction
%
% model parameters:
%   R_type - which relations are used in the decomposition; modify this
%   setting to run MetaFac digg prediction ('R14') / comment prediction ('R15'),
%   or MWA digg prediction ('FullD') / comment prediction ('FullC')
% other configuration:
%   other model / data parameters not specified in this file can be
%   configured in get_digg_option.m

if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);

R_types = {'R2','R23','R123','R234','R1234','R2345','R12345','R23456','R123456', ...
    'R13','R14','R15','R134','R135','R145','R1345', ...
 'FullP','FullD','FullC','FullPCD', ...
 'R4789','R125','R5789','R478','R124', ...
 'R1235','R1245'};
% R=11 - R14 basic digg prediction
% R=12 - R15 basic comment prediction
% R=18,19 - baseline multi-way aspect for digg / comment prediction

filename=sprintf('%sblog_%s_v%d-%d_%s%d%s%d',outfilepath,'rntf_foldin_predict_diggs_story',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename '.txt'];
savelog(logfile,logfile);
opts.foldin_sz = 1;

% for R = [ 11 12 22 14 26 15 27 16 7 9]
for R = 11
    opts.R_type=R_types{R};    
    for K = 4 %[4:4:20]
        opts.K=K;
     for pa = 0.2 %0:0.2:1
        opts.pa = pa;
        if R==18 || R==19,opts.pa=0;end
        s = sprintf('predict in R_type=%s K=%d pa=%f...',opts.R_type,opts.K,opts.pa);
        savelog(logfile,s);
        if R<=16 || R>20
            run_diggs_rntf(opts);
        else
            run_diggs_full_rntf(opts);
        end
        foldin_diggs_story2(opts);     
        predict_diggs_story(opts);
        if R==18 || R==19,break;end
     end
    end
    %}
end