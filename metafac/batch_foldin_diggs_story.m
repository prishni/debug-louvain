function batch_foldin_diggs_story()
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);

R_types = {'R2','R23','R123','R234','R1234','R2345','R12345','R23456','R123456', ...
    'R13','R14','R15','R134','R135','R145','R1345', ...
 'FullP','FullD','FullC','FullDCP', ...
 'R4789'};

filename=sprintf('%sblog_%s_v%d-%d_%s%d%s%d',outfilepath,'foldin_diggs_story',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename '.txt'];
savelog(logfile,logfile);
opts.foldin_sz = 1;

for R=21 %[11 16] %[11 13 15 16]
    opts.R_type=R_types{R};
    for K = 4:4:20
%     for K = [2:2:20 30:10:50]
%     for K = [4 12] 
    opts.K=K;
%     for pa = 0.1:0.2:0.9
        opts.pa = pa;
        s = sprintf('folding in R_type=%s K=%d pa=%f ...',opts.R_type,opts.K,opts.pa);
        savelog(logfile,s);
        foldin_diggs_story2(opts);
%     end
    end
end