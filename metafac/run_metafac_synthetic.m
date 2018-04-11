for i = 18:20
    for p = [8]
        for a = [4,6,8,10]
            for mu = [40]
                for d = [4]
                    f00 = ['C:\Users\Home\Desktop\Anchit2\MetaFac\' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(d) '_00'];
                    f11 = ['C:\Users\Home\Desktop\Anchit2\MetaFac\' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(d) '_11'];
                    f01 = ['C:\Users\Home\Desktop\Anchit2\MetaFac\' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(d) '_01'];
                    if exist(f00, 'file') ~= 2
                        disp({p, a, mu, d})
                        disp('LOL')
                        continue
                    end
                    f_00 = fopen(f00, 'r');
                    f_11 = fopen(f11, 'r');
                    f_01 = fopen(f01, 'r');
                    fd_00 = fscanf(f_00, '%d %d', [2 Inf])+1;
                    fd_11 = fscanf(f_11, '%d %d', [2 Inf])+1;
                    fd_01 = fscanf(f_01, '%d %d', [2 Inf])+1;
                    fd_00 = sort(fd_00', 2);
                    fs_00 = sptensor(fd_00, 1);
                    if sum(size(fs_00) ~= [100 100]) ~= 0
                        fs_00(100, 100) = 0
                    end

                    fd_11 = sort(fd_11', 2);
                    z = ones(size(fd_11))*100;
                    fd_11 = fd_11 - z;
                    fs_11 = sptensor(fd_11, 1);
                    if sum(size(fs_11) ~= [100 100]) ~= 0
                        fs_11(100, 100) = 0
                    end

                    fd_01 = sort(fd_01', 2);
                    z = zeros(size(fd_01));
                    z(:, 2) = 100;
                    fd_01 = fd_01 - z;
                    fs_01 = sptensor(fd_01, 1);
                    if sum(size(fs_01) ~= [100 100]) ~= 0
                        fs_01(100, 100) = 0
                    end

                    RG = {{fs_00, [1 1], 1}, {fs_11, [2 2], 1}, {fs_01, [1 2], 1}};
                    Vdims = [100 100];
                    for K = [5:3:30] 
                        fK = ['C:\Users\Home\Desktop\Anchit2\MetaFacRes\' num2str(i) '_' num2str(p) '_' num2str(a) '_' num2str(mu) '_' num2str(d) num2str(K) '.mat'];
                        [S, U, iters, cvtime, ll] = rntf(RG, Vdims, K, {});
                        save(fK, 'S', 'U', 'iters', 'cvtime', 'll')
                    end
                end
            end
        end
    end
end
