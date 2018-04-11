uc = 244
lc = 1627
for K = [5:3:30]
    fK = ['C:\Users\Home\Desktop\Anchit2\MetaFacRes4\' num2str(uc) '_' num2str(K) '.mat'];
    load(fK);
    communities = zeros(K,uc+lc);
    for m = [1:uc]
        [dummy, mid] = max(U{1}(m,:) .* S');
        communities(mid, m) = 1;
    end
    for m = [1:lc]
        [dummy, mid] = max(U{2}(m,:) .* S');
        communities(mid, m+uc) = 1;
    end
    fS = ['C:\Users\Home\Desktop\Anchit2\MetaFacCom4\' num2str(uc) '_' num2str(K)];
    fileID = fopen(fS, 'w');
    for c = [1:K]
        for d = [1:uc+lc]
            if communities(c, d) == 1
                fprintf(fileID, '%d ', d);
            end
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end