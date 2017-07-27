clear all
close all
clc

%Problem params
n_img = 3;
n_sz = 3;
n_filter = 3;
N = n_img*n_sz*n_filter;
n_rep = 10;

%Function handles
f1 = @(I,n) conv2(I,fspecial('average',n),'same');
f2 = @(I,n) conv2(I,fspecial('gaussian',n),'same');
f3 = @(I,n) wiener2(I,[n n]);
funcs = {f1,f2,f3};

%Window sizes
sz = [3 5 7];

%Load images and convert to grayscale
imgs_orig = cell(n_img,1);
imgs = cell(n_img,1);
%original_psnr = zeros(n_img,1);
d = dir('*.jpg');
for i=1:n_img
    aux = imread(d(i).name);
    if(size(aux,3) > 1)
        aux = rgb2gray(aux);
    end
    imgs_orig{i} = mat2gray(aux);
    %imgs{i} = double(aux);
    %imgs{i} = imnoise(aux,'gaussian'); %0 mean, 0.01 var
    %original_psnr(i) = psnr(uint8(imgs{i}),uint8(aux));
    %imgs{i} = double(imgs{i});
end

%Result storage
original_psnr = zeros(n_img,1);
res_psnr = zeros(n_img,n_sz,n_filter,n_rep);
res_t = zeros(n_img,n_sz,n_filter,n_rep);
order = zeros(N*n_rep,1);

%Generate all possible combinations
% [B,C,A] = meshgrid(1:n_img,1:n_sz,1:n_filter);
% combs_aux = [reshape(A,[1 numel(A)]) ; reshape(B,[1 numel(B)]) ; reshape(C,[1 numel(C)])]';
combs_aux = fullfact([n_img n_filter n_sz]);

for i=1:n_rep
    combs = combs_aux; %Populate the combination array again
    for j=1:n_img
        imgs{j} = imnoise(imgs_orig{j},'gaussian'); %Add noise to img
        original_psnr(j) = psnr(imgs{j},imgs_orig{j});
    end
    for j=1:N
        idx = randi(size(combs,1),1);
        order((i-1)*N + j) = idx;
        p = combs(idx,:); %[img #, filter size, filter type]
        combs(idx,:) = [];
        aux = funcs{p(3)};
        %img = imnoise(imgs{p(1)},'gaussian');
        %original_psnr = psnr(img,imgs{p(1)});
        %img = imgs{p(1)};
        tic;
        filtered = aux(imgs{p(1)},sz(p(2)));
        res_t(p(1),p(2),p(3),i) = toc;
        res_psnr(p(1),p(2),p(3),i) = psnr(filtered,imgs_orig{p(1)});
    end
    res_psnr(:,:,:,i) = res_psnr(:,:,:,i)./repmat(original_psnr,[1 3 3]);
end

res_t = res_t.*1000; %Convert to ms
res_vector = reshape(res_psnr,[1 numel(res_psnr)]);
times_vector = reshape(res_t,[1 numel(res_t)]);
order_vector = repmat(combs_aux',[1 n_rep]);
% p = anovan(res_vector,{order_vector(1,:),order_vector(2,:),order_vector(3,:)},'model','interaction','varnames',{'Img','Size','Type'});
