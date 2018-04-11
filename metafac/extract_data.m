for i = 11:20
    for p = [8]
        for a = [4,6,8,10]
            for mu = [40]
                for den = [4]
                    for K = [5:3:30]
                        fK = ['./MetaFacRes/' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(den) num2str(K) '.mat'];
                        if exist(fK, 'file') ~= 2
                            disp({p, a, mu, den})
                            disp('LOL')
                            continue
                        end
                        load(fK);
                        communities = zeros(K,200);
                        for u = [1 2]
                            for m = [1:100]
                                [dummy, mid] = max(U{u}(m,:) .* S');
                                if u == 1
                                    communities(mid, m) = 1;
                                else
                                    communities(mid, m+100) = 1;
                                end
                            end
                        end
                        fS = ['./MetaFacCom/' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(den) num2str(K)];
                        fileID = fopen(fS, 'w');
                        for c = [1:K]
                            for d = [1:200]
                                if communities(c, d) == 1
                                    fprintf(fileID, '%d ', d);
                                end
                            end
                            fprintf(fileID, '\n');
                        end
                        fclose(fileID);
                    end
                end
            end
        end
    end
end