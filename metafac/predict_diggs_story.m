function predict_diggs_story(opts)
% compute p(s|z)p(u|z) and output ordered list

% model params
if nargin == 0,opts = struct();end
[K,pa,TMIN,TMAX,ndt,R_type,infilepath0,infilepath,outfilepath,data_ver,res_ver,stream,stream_s,topN,topP,foldin_sz] = get_digg_option(opts);
discount = 1./log2(1+(1:topP));
perf = [];
for ti = 2:TMAX
    fprintf('\nt=%d',ti);
    file_prefix=sprintf('%sfoldin_story_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
    filename = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 'f' num2str(foldin_sz) 't' num2str(ti) '.mat'];
    load(filename,'U_s','sid'); 
    nNS = length(sid);

    filename=sprintf('%starget_%s%d_v%d.mat',infilepath,stream_s,ti,data_ver);
    load(filename,'P','D','C','R','nU','nS');
    switch R_type
        case {'FullD','R14','R4789','R478','R124','R1234'}
            A=D;
        case {'FullC','R15','R125','R5789','R135','R145','R1235','R1245','R1345','R12345','R123456'}
            A=C;
        case {'FullDCP'}
            A=D+C+P; %+R
    end
    A = A(:,sid);
    uid = find(sum(A,2)); % find users who act on give story list
    nNU = length(uid);
    A = A(uid,:)>0; %nnz(A)
    ps = sum(A,1); ps = ps./ sum(ps);

    file_prefix=sprintf('%sres_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
    decom_res_file = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 't' num2str(ti-1) '.mat'];
    if exist(decom_res_file,'file')>0
        load(decom_res_file,'S','U');  
        fprintf('\nload %s',decom_res_file);
    else
        error('\nfile not found:%s',decom_res_file)
    end    

    switch R_type
        case {'FullD','FullC','FullP','FullDCP'}
        U_u=U{1};
        otherwise
        faid = get_diggs_facet_id(R_type,'user');
        U_u = U{faid}(uid,:);
    end
    % customer prediction
%     pred_i = sparse(nNU,topP);
%     pred_v = sparse(nNU,topP);
%     pred_t = sparse(nNU,topP);
    hitc = 0; 
    for i=1:nNU
        pr = sparse(1,nNS);
        for z =1:K
            if foldin_sz
            pr = pr + U_u(i,z).*S(z).*U_s(:,z)';
            else
            pr = pr + U_u(i,z).*U_s(:,z)';
            end
        end % z
        us = pr .* ps; % multiply by p(story)
        [v ind]=sort(-us);
        ind = ind(1:topP);
%         v = -v(1:topP);
%         pred_i(i,:) = ind;
%         pred_v(i,:) = v;
%         pred_t(i,:) = A(i,ind);
        pred_idx = ind; target_idx = find(A(i,:));
        [hit pi targeti]=intersect(pred_idx, target_idx);

        prec(i) = length(hit)/topP;
        recl(i) = length(hit)/length(target_idx);
        cg(i) = length(hit);
        bestn = min([topP length(target_idx)]);
        maxcg(i) = bestn;
        dcg(i) = sum(discount(pi));
        best = 1:bestn;
        maxdcg(i) = sum(discount(best));
        nnru(i) = length(target_idx);
        hitcnt(i)=length(hit)>0;
%         if mod(i,500) == 0,fprintf('\nprogress: %d ',i);end
    end
    nnu=nNU; nnr=nNS;
    nnru = sum(nnru)/nnu;
    prec = sum(prec)/nnu;
    recl = sum(recl)/nnu;
    ncg = sum(cg)/sum(maxcg);
    ndcg = sum(dcg)/sum(maxdcg);
    hitcnt = sum(hitcnt)/nnu;
    fprintf('\nperformance: %.3f %.3f %.3f %.3f %.3f (nnu=%d,nnr=%d)',prec,recl,ncg,ndcg,hitcnt,nnu,nnr);

    % global prediction
%     [ii jj v] = find(pred_v);
%     [v ind] = sort(-v); length(ind)
%     li = sub2ind([nNU,topP],ii(ind),jj(ind));
%     size(pred_t(li)),input('.')
%     pred_g = [ii(ind), jj(ind), -v, pred_t(li)];
    
    perf = [perf; hitcnt ndcg ncg prec recl nnu nnr];
    file_prefix=sprintf('%spredres_story_v%d-%d_%s%d%s%d',outfilepath,data_ver,res_ver,stream_s,TMIN,stream_s,TMAX);
    filename = [file_prefix 'K' num2str(K) 'pa' num2str(pa) R_type 'f' num2str(foldin_sz) '.mat'];
    save(filename,'perf');
    
end % ti