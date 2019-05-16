close all; clearvars -except set_up
%% Set Up
addpath(genpath('functions'));
FFW = set_up.FFW;
FB = set_up.FB;
% Take a look at the pdf -> 106 
ID_FFW = EI3NodeMotif.weimat2ID(FFW);
% Take a look at the pdf -> 114
ID_FB = EI3NodeMotif.weimat2ID(FB);
%% Test Matrix
TM = set_up.TM; %test matrix
N = length(TM); %number of vertices
for i=1:N
    for j=1:N
        TM(i,j)=randsample(-1:1,1);
        if i==j
            TM(i,j) = 0;
        end
    end
end
%% Random Subgraph Sampling
C = nchoosek(N,3); %factorial(N)/(factorial(3)*factorial(N-3)); number of all combinations
comb3node = nchoosek(1:N,3); %all possible combinations
percent = set_up.percent; %percentage of sampling from combinations
pc = percent*C;
%% Iterations
iter = set_up.iter;
fprintf('Iterations...')
tic;
for q = 1:iter
    ID_TM = zeros(size(pc,1),1);
    c = cell(pc,1); %initiate cell array for select combinations
    sM = cell(pc,1); %sampled matrices
    sV = cell(pc,1); %matrix indexes (upper symmetric, lower symmetric) w/o diagonal
    for i = 1:pc
        c{i,1} = randsample(N,3); %=subsampled comb3node (in Tuan's code)
        c{i,1} = sort(c{i,1});
        sM{i,1} = TM([c{i,1}],[c{i,1}]); %=g_i (in Tuan's code)
        sV{i,1} = sM{i,1}([4,7,8 2,3,6]);
        ID_TM(i,1) = EI3NodeMotif.weimat2ID(sM{i,1});
    end
    %% Assigning Unique IDs
    unq_ids = unique(ID_TM); %sort unique IDs
    counts = arrayfun(@(x) sum(ID_TM == x),unq_ids,'uni',1); %corresponding counts
    check=sum(counts)==pc; %if 1, no errors
    if check ~=1
        msg = 'error: number not add up';
    end
    ffwlook4 = 106; %looking for FFW
    fblook4 = 114; %looking for FB
    %% Finding Subgraphs with Specific IDs
    j=1; k=1;
    for i = 1:size(ID_TM,1)
        if ID_TM(i,1) == ffwlook4
            w(j,1) = i;
            j = j+1;
        elseif ID_TM(i,1) == fblook4
            b(k,1) = i;
            k = k+1;
        end
    end

    if j >1
        ffw = cell(length(w),1);
        for i = 1:length(w)
            ffw{i,1} = c{w(i,1),1};
        end 
    end
    if k > 1
        fb = cell(length(b),1);
        for i = 1:length(b)
            fb{i,1} = c{b(i,1),1};
        end
    end
    %% Ratio of FFW and FB
    nffw = size(ffw,1);
    nfb = size(fb,1);
    ratio(q,1) = nffw/nfb;
end
ratio_m = mean(ratio);
ratio_std = sqrt(var(ratio));

iter_time = toc;
fprintf('took %.3f seconds. \n', iter_time)