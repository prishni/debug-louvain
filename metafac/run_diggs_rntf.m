function run_diggs_rntf(opts)
% run hypergraph factorization on Diggs data


% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s] = get_digg_option(opts);
filename=sprintf('%slog_%s_v%d-%d_%s%d%s%d',outfilepath,'run_diggs_rntf',data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
logfile = [filename 'K' num2str(K) R_type 'pa' num2str(pa) '.txt'];
lfile = fopen(logfile, 'w');
fprintf(lfile,logfile);
file_prefix=sprintf('%sres_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);

%load static data tensors
filename=sprintf('%ssUC_v%d.mat',infilepath,data_ver); 
load(filename); sUC=sptensor(sUC);% get sUC (user,contact)
% filename=sprintf('%ssW_v%d.mat',infilepath,data_ver);    
% load(filename); sW=sptensor(sW); % get sW (story,keyword,topic)

tstart=TMIN; tmax=TMAX;
for ti = 1:tmax
    if ti<tstart,continue;end  
    % load data tensors
    filename=sprintf('%ssP%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename); sP=sptensor(sP); % get sP (user,story) via submit
    filename=sprintf('%ssD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename); sD=sptensor(sD);% get sD (user,story) via digg
    filename=sprintf('%ssC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename); sC=sptensor(sC);% get sC (user,story,comment) 
    filename=sprintf('%ssR%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename); sR=sptensor(sR);% get sR (user,story,comment) 
    filename=sprintf('%ssW%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename); sW=sptensor(sW);% get sW (story,keyword,topic) 
%     for tt=1:ti-1
%         filename=sprintf('%ssW%s%d_v%d.mat',infilepath,stream_s,tt,data_ver);
%         load(filename); sW=sW+sptensor(sW);% get sW (story,keyword,topic) 
%     end

    switch R_type
        case {'R4789','R5789','R478'}
        filename=sprintf('%ssSW%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sSW=sptensor(sSW);% get sSW (story,keyword) 
        filename=sprintf('%ssSJ%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sSJ=sptensor(sSJ);% get sSJ (story,topic) 
        filename=sprintf('%ssWJ%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sWJ=sptensor(sWJ);% get sWJ (keyword,topic) 
        sW=[];
        otherwise
            sSW=[];sSJ=[];sWJ=[];
    end    
    all_dims = [nU,nS,nC,nW,nJ];
%     all_tensors = {sW,sUC,sP,sD,sC,sR};
    fprintf('\nt=%d',ti);
    if ti>tstart,prior = {S,U,pa};
    else prior={}; end
    
    decom_res_file = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 't' num2str(ti) '.mat'];
    if exist(decom_res_file,'file')>0
        load(decom_res_file,'S','U');  
        fprintf('\n%s exist',decom_res_file);
        fprintf(lfile,sprintf('%s exist',decom_res_file));
    else
%         [Vdims RG]=mk_hg_tensor(R_type,all_dims,sW,sUC,sP,sD,sC,sR);
        [Vdims RG]=mk_hg_tensor(R_type,all_dims,sW,sUC,sP,sD,sC,sR,sSW,sSJ,sWJ);
        [S,U,iters,cvtime,ll]=rntf(RG,Vdims,K,prior);
        save(decom_res_file,'S','U','iters','cvtime','ll');  
        fprintf('\nsave %s',decom_res_file);
        fprintf(lfile, sprintf('save %s',decom_res_file));    
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Vdims RG] = mk_hg_tensor(R_type,all_dims,sW,sUC,sP,sD,sC,sR,sSW,sSJ,sWJ)
    switch(R_type)
        case {'R2'} % user-contact (contact)
            Vdims = [all_dims(1)];
            RG = {{sUC,[1 1],1}};
        case {'R23'} % contact,user-story (submit)
            Vdims = [all_dims(1:2)];
            RG = {{sUC,[1 1],1},{sP,[1 2],1}};
        case {'R123'} % story-keyword-topic (topic),contact,submit
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sW,[2 3 4],1},{sUC,[1 1],1},{sP,[1 2],1}};
        case {'R234'} % contact,submit,user-story (digg)
            Vdims = [all_dims(1:2)];
            RG = {{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1}};
        case {'R1234'} % topic, contact,submit,digg
            Vdims = [all_dims([1:2 4:5]) ];
            RG = {{sW,[2 3 4],1},{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1}};
        case {'R2345'} % contact,submit,digg,user-story-comment (comment)
            Vdims = [all_dims(1:3)];
            RG = {{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1},{sC,[1 2 3],1}};
        case {'R12345'} % topic,contact,submit,digg,user-story-comment (comment)
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1},{sC,[1 2 3],1}};
        case {'R23456'} % contact,submit,digg,comment,user-comment-reply (reply)
            Vdims = [all_dims(1:3)];
            RG = {{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1},{sC,[1 2 3],1},{sR,[1,3],1}};
        case {'R123456'} % topic,contact,submit,digg,comment,user-comment-reply (reply)
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sUC,[1 1],1},{sP,[1 2],1},{sD,[1 2],1},{sC,[1 2 3],1},{sR,[1,3],1}};
        case {'R13'} % topic,submit
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sW,[2 3 4],1},{sP,[1 2],1}};
        case {'R14'} % topic,digg
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sW,[2 3 4],1},{sD,[1 2],1}};
        case {'R15'} % topic,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sC,[1 2 3],1}};
        case {'R125'} % topic,comment,contact
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sC,[1 2 3],1},{sUC,[1 1],1}};
        case {'R134'} % topic,submit,digg
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sW,[2 3 4],1},{sP,[1 2],1},{sD,[1 2],1}};
        case {'R135'} % topic,submit,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sP,[1 2],1},{sC,[1 2 3],1}};
        case {'R1235'} % topic,contact,submit,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sP,[1 2],1},{sC,[1 2 3],1},{sUC,[1 1],1}};
        case {'R145'} % topic,digg,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sD,[1 2],1},{sC,[1 2 3],1},{sUC,[1 1],1}};
        case {'R1245'} % topic,digg,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sD,[1 2],1},{sC,[1 2 3],1}};
        case {'R1345'} % topic,submit,digg,comment
            Vdims = [all_dims(1:5)];
            RG = {{sW,[2 4 5],1},{sP,[1 2],1},{sD,[1 2],1},{sC,[1 2 3],1}};
        case {'R4789'} % topic (story-keyword,story-topic),digg
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sD,[1 2],1},{sSW,[2 3],1},{sSJ,[2 4],1},{sWJ,[3 4],1}};            
        case {'R478'} % topic (story-keyword,story-topic),digg
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sD,[1 2],1},{sSW,[2 3],1},{sSJ,[2 4],1}};            
        case {'R124'} % story-keyword-topic (topic),contact,submit
            Vdims = [all_dims([1:2 4:5])];
            RG = {{sW,[2 3 4],1},{sUC,[1 1],1},{sD,[1 2],1}};            
        case {'R5789'} % topic,comment
            Vdims = [all_dims(1:5)];
            RG = {{sC,[1 2 3],1},{sSW,[2 4],1},{sSJ,[2 5],1},{sWJ,[4 5],1}};
    end