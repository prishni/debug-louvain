f00 = ['C:\Users\Home\Desktop\Anchit2\MetaFac4\244_1627_2\00'];
f11 = ['C:\Users\Home\Desktop\Anchit2\MetaFac4\244_1627_2\11'];
f01 = ['C:\Users\Home\Desktop\Anchit2\MetaFac4\244_1627_2\01'];
if exist(f00, 'file') ~= 2
    disp({p, a, mu, d})
    disp('LOL')
    continue
end

uc = 244
lc = 1627
f_00 = fopen(f00, 'r');
f_11 = fopen(f11, 'r');
f_01 = fopen(f01, 'r');
fd_00 = fscanf(f_00, '%d %d', [2 Inf])+1;
fd_11 = fscanf(f_11, '%d %d', [2 Inf])+1;
fd_01 = fscanf(f_01, '%d %d', [2 Inf])+1;
fd_00 = sort(fd_00', 2);
fs_00 = sptensor(fd_00, 1);
if sum(size(fs_00) ~= [uc uc]) ~= 0
    fs_00(uc, uc) = 0
end

fd_11 = sort(fd_11', 2);
z = ones(size(fd_11))*uc;
fd_11 = fd_11 - z;
fs_11 = sptensor(fd_11, 1);
if sum(size(fs_11) ~= [lc lc]) ~= 0
    fs_11(lc, lc) = 0
end

fd_01 = sort(fd_01', 2);
z = zeros(size(fd_01));
z(:, 2) = uc;
fd_01 = fd_01 - z;
fs_01 = sptensor(fd_01, 1);
if sum(size(fs_01) ~= [uc, lc]) ~= 0
    fs_01(uc, lc) = 0
end

RG = {{fs_00, [1 1], 1}, {fs_11, [2 2], 1}, {fs_01, [1 2], 1}};
Vdims = [uc lc];
for K = [5:3:30] 
    fK = ['C:\Users\Home\Desktop\Anchit2\MetaFacRes4\244_' num2str(K) '.mat'];
    [S, U, iters, cvtime, ll] = rntf(RG, Vdims, K, {});
    save(fK, 'S', 'U', 'iters', 'cvtime', 'll')
end
-