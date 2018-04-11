function [U_s sid U_u S] = foldin_diggs_story(opts,ti)
% compute p(s|z) by folding-in

% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);

outfile_prefix=sprintf('%sfoldin_story_PARAFAC_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);

% for ti = 2:TMAX
    fprintf('\nt=%d',ti);
    switch R_type
        case {'FullD'}
        filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = sFullD; clear sFullD;
        case {'FullC'}
        filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullC=sptensor(sFullC);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = sFullC; clear sFullC;
        case {'FullDCP'}
        filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = sFullD; clear sFullD;
        filename=sprintf('%ssFullC%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = nzS+sFullD; clear sFullD;
        filename=sprintf('%ssFullP%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullP=sptensor(sFullP);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = nzS+sFullP; clear sFullP;
    end
    sid = collapse(nzS, [1 3 4]);
    sid = find(sid); 
    nNS=length(sid);
    sid2i =sparse(1,nS); % inverse mapping from original sid to 1:nNS
    for i=1:nNS,sid2i(sid(i))=i;end

    file_prefix=sprintf('%sresPARAFAC_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
    decom_res_file = [file_prefix 'K' num2str(K) R_type 't' num2str(ti-1) '.mat'] ;
    if exist(decom_res_file,'file')>0
        load(decom_res_file,'S','U');  
        fprintf('\nload %s',decom_res_file);
    end    

    U_u = U{1};
    U_w = U{3};
    U_j = U{4};
    [subs vals] = find(nzS); 
    uu = subs(:,1);
    ss = subs(:,2);
    ss=sid2i(ss);
    ww = subs(:,3);
    jj = subs(:,4);
    U_s = sparse(nNS,K);
    
    for z=1:K
        U = U_u(uu,z).*U_w(ww,z).*U_j(jj,z)./S(z);
        U = collapse(sptensor(subs,U), [1 3 4]);
        U_s(:,z) = U(sid);
    end
% end % ti

