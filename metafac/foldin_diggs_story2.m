function foldin_diggs_story(opts)
% compute p(s|z) by folding-in

% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);

outfile_prefix=sprintf('%sfoldin_story_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);

for ti = 2:TMAX
    fprintf('\nt=%d',ti);
    switch R_type
        case {'FullD','R14','R4789','R478','R124','R1234'}
        filename=sprintf('%ssFullD%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
        load(filename); sFullD=sptensor(sFullD);% get sD (user,story) via digg
        % (user,story,keyword,topic)
        nzS = sFullD; clear sFullD;
        case {'FullC','R15','R125','R5789','R135','R1235','R145','R1245','R1345','R12345','R123456'}
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

    file_prefix=sprintf('%sres_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
    decom_res_file = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 't' num2str(ti-1) '.mat'];
    if exist(decom_res_file,'file')>0
        load(decom_res_file,'S','U');  
        fprintf('\nload %s',decom_res_file);
    else
        error('\nfile not found:%s',decom_res_file)
    end    

    switch R_type
        case {'FullD','FullC','FullDCP'}
        U_u = U{1}; U_j = U{4}; U_w = U{3};
        otherwise
        faid = get_diggs_facet_id(R_type,'user');
        U_u = U{faid};
        faid = get_diggs_facet_id(R_type,'topic');
        U_j = U{faid};
        faid = get_diggs_facet_id(R_type,'keyword');
        U_w = U{faid};
    end
    % initialize U_s
    nK = size(U_u,2); % no. of factor
    U_s = rand(nNS, nK)+0.5;
    if foldin_sz
        nrzU = sum(U_s); U_s = U_s ./ repmat(nrzU, [nNS,1]); % p(s|z)
    else
        nrzU = sum(U_s,2); U_s = U_s ./ repmat(nrzU, [1,nK]); % p(z|s)
    end

    [subs vals] = find(nzS); 
    uu = subs(:,1);
    ss = subs(:,2);
    ss=sid2i(ss);
    ww = subs(:,3);
    jj = subs(:,4);

    maxiter=500; minper=1e-6;
    logl = zeros(1,maxiter);
    tic;
    for iter=1:maxiter
        zzS = sparse(length(vals),1);
        % M-step
        for z=1:nK
            if foldin_sz
                % multiply S(z) might tend to overfitting, not obvious... 
                % use this in archival-pc training (20090127)
            zS{z} = U_s(ss,z).*U_u(uu,z).*U_w(ww,z).*U_j(jj,z).*S(z) + (1e-20);
            else
            zS{z} = U_s(ss,z).*U_u(uu,z).*U_w(ww,z).*U_j(jj,z) + (1e-20);
            end
            zzS = zzS + zS{z};
        end % z
        zzS(find(zzS==0))=1;
        for z = 1:nK
            zS{z} = zS{z} ./ zzS;
            zS{z} = sptensor(subs,zS{z});
        end %z

        % E-step
        U_s = zeros(nNS,nK);
        for z = 1:nK
            v = double(collapse(zS{z},[1 3 4]));  
            U_s(:,z) = v(sid);
        end %z
        if foldin_sz
            nrzU = sum(U_s); U_s = U_s ./ repmat(nrzU, [nNS,1]); 
        else
            nrzU = sum(U_s,2); nrzU(find(nrzU==0))=1; U_s = U_s ./ repmat(nrzU, [1,nK]); 
        end

        ll=sum(vals.*log(zzS));
        logl(iter) = ll;
        if iter>1 && abs((logl(iter)-logl(iter-1))/logl(iter-1)) < minper,
            break;
        end
    end % iter
    cvtime=toc;
    foldin_file = [outfile_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 'f' num2str(foldin_sz) 't' num2str(ti) '.mat'];
    save(foldin_file,'U_s','sid','sid2i','iter','cvtime');
end % ti

