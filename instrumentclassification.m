function J_spat=instrumentclassification(music)
results_dir  = 'output_data/';
output_prefix='output';

stft_win_len = 2048;
qerb_nbin = 350;       % number of frequency bins for qerb

J_spat = 4;       % number of spatial components
J_spec_spat = 4;  % number of spectral components per one spatial component

EM_iter_num = 50;

% add '/' character at the end of directory names
% Input time-frequency representation
fprintf('Input time-frequency representation\n');
[x, fs]=audioread(music);
x = x.';
mix_nsamp = size(x,2);

% compute time-frequency representation
Cx = comp_transf_Cx(x, 'qerb', stft_win_len, fs, qerb_nbin);

[F, N, ~, M] = size(Cx);

% construction of wide harmonic blobs
hcomp=1;
wcomp=0;
hwidth=44/3;
wwidth=44/3;
eslope=6;
firstnote=21;
lastnote=108;

[V_harm,B_harm]=spectral_blobs_erb(qerb_nbin,fs,hcomp,wcomp,firstnote,lastnote,hwidth,wwidth,eslope);

wcomp = length(find(max(B_harm, [], 2) == 0));

V_harm = V_harm(:, 1:end-wcomp);
B_harm = B_harm(1:end-wcomp, :);

% load drums spectral shapes
load('spectral_patterns/drums_spec_shepes_clust_GSMM.mat', 'W_clust_GSMM');

W_drums = W_clust_GSMM;

% load bass spectral shapes
load('spectral_patterns/bass_spec_shepes_clust_GSMM.mat', 'W_clust_GSMM');

% W_bass = W_clust_GSMM / 8;

mix_str.Cx = Cx;
mix_str.transf = 'qerb';
mix_str.Noise_PSD = zeros(F, 1);
mix_str.spat_comps = cell(1,J_spat);
mix_str.spec_comps = cell(1,J_spat * J_spec_spat);
mix_str.fs = fs;
mix_str.wlen = stft_win_len;

spec_ind = 1;
for spat_ind = 1:J_spat
    % initialize spatial component
    spat_comp = mix_str.spat_comps{spat_ind};
    
    rank = 2;
    spat_comp.time_dep   = 'indep';
    spat_comp.mix_type   = 'conv';
    spat_comp.frdm_prior = 'free';

    % initialize parameters of filters
    if strcmp(spat_comp.mix_type, 'inst')
        params = ones(M,rank,1);
    else
        params = ones(M,rank,F);
    end;

    params = 0.5 * (1.8 * abs(randn(size(params))) + 0.2 * ones(size(params)));
    
    if strcmp(spat_comp.mix_type, 'conv')
        params = params .* sign(randn(size(params)) + sqrt(-1)*randn(size(params)));
        for f = 1:F
            params(:,:,f) = params(:,:,1);
        end;
    end;    

    % adjust total parameters enegry
    spat_comp.params = params * sqrt(mean(abs(mix_str.Cx(:))));

    mix_str.spat_comps{spat_ind} = spat_comp;
    
    % initialize spectral components
    for j = 1:J_spec_spat
        spec_comp = mix_str.spec_comps{spec_ind};
        
        spec_comp.factors = cell(1, 1);

        spec_comp.spat_comp_ind = spat_ind;

        % initialize one factor
        factor = spec_comp.factors{1};
        
        switch j
           case 1
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];

                factor.TW_constr = 'NMF';
                
                factor.FB = V_harm;
                factor.FW = B_harm;

                K = size(factor.FW, 2);
           case 2
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];

                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);

                factor.FW = diag(ones(1, K));
           case 3
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
                
           case 4
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
           case 5
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
           case 6
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
           case 7
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
           case 8
                factor.FB_frdm_prior = 'fixed';
                factor.FW_frdm_prior = 'free';
                factor.TW_frdm_prior = 'free';
                factor.TB_frdm_prior = [];
                
                factor.TW_constr = 'NMF';
                
                factor.FB = W_drums;
                
                K = size(factor.FB, 2);
                
                factor.FW = diag(ones(1, K));
               
                
           otherwise
                disp('Variable j must be smaller than 3.')
        end;
        
        factor.TW = 0.5 * (1.5 * abs(randn(K, N)) + 0.5 * ones(K, N));
        factor.TB = [];
        
        spec_comp.factors{1} = factor;
        
        mix_str.spec_comps{spec_ind} = spec_comp;
    
        spec_ind = spec_ind + 1;
    end;
end;

sig_pow = zeros(F, N);
for m = 1:M
    sig_pow = sig_pow + mix_str.Cx(:,:,m,m) / M;
end;
Ann_PSD_beg = mean(sig_pow, 2) / 10;

% run EM algorithm for parameters estimation
mix_str = estim_param_a_post_model(mix_str, EM_iter_num, 'ann', Ann_PSD_beg);

% determine the predominant harmonic source
prev_mean_TW = 0;
pred_ind = 1;

spec_ind = 1;
for spat_ind = 1:J_spat
    for j = 1:J_spec_spat
        if j == 1
            spec_comp_factor1 = mix_str.spec_comps{spec_ind}.factors{1};

            mean_TW = mean(spec_comp_factor1.TW(:));

            if prev_mean_TW < mean_TW
                prev_mean_TW = mean_TW;
                pred_ind = spec_ind;
            end;
        end;
        
        spec_ind = spec_ind + 1;
    end;
end;

% fill in the separated components indices
sep_cmp_inds = cell(1, J_spec_spat);
spec_ind = 1;
for spat_ind = 1:J_spat
    for j = 1:J_spec_spat
        if j == 1
            if spec_ind == pred_ind
                sep_cmp_inds{1} = [sep_cmp_inds{1}, spec_ind];
            else
                sep_cmp_inds{4} = [sep_cmp_inds{4}, spec_ind];
            end;
        else
            sep_cmp_inds{j} = [sep_cmp_inds{j}, spec_ind];
        end;
        
        spec_ind = spec_ind + 1;
    end;
end;

i_est = separate_spec_comps(x, mix_str, sep_cmp_inds);

for i=1:J_spat
    audiowrite([results_dir output_prefix strcat(num2str(i),'.wav')],reshape(i_est(i,:,:),mix_nsamp,2),fs);
end