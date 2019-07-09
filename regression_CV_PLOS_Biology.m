clear all
close all
cd /Users/administrator/Documents/Postdoc/Phenotyping_personality2/
text1 = textread(['/Users/administrator/Documents/Postdoc/Phenotyping_personality2/Matlab/List/Blind_Subjects_V2.txt'],'%s'); % naming and loading the subject list file
Pers_V2{1} = importdata('/Users/administrator/Documents/Postdoc/Phenotyping_personality2/Matlab/Data/pain_trait.mat');
Pers_V2{2} = importdata('/Users/administrator/Documents/Postdoc/Phenotyping_personality2/Matlab/Data/emote_trait.mat');

load('/Users/administrator/Documents/Postdoc/Phenotyping_personality2/Matlab/Data/FD.mat');
mFD = mean(FD,2);

folds_d{1} = [21:61];
folds_d{2} = [1:20 41:61];
folds_d{3} = [1:40];

folds_t{1} = [1:20];
folds_t{2} = [21:40];
folds_t{3} = [41:61];

ns = length(text1); % determine the number of subjects

%% organize connectome
for l = 1:size(Pers_V2,2)
for k = 1:size(folds_d,2)
for n =1:ns;
adj_mat = importdata(['/Users/administrator/Documents/Data/Placebo_1/matrices/visit2/total/filt_MWVG_' text1{n} '_adj_matrix_limbic.txt']);
r=adj_mat(:);
z=.5.*log((1+r)./(1-r));
connectome(n,:) = z;
end

discovery = connectome(folds_d{k},:)';
perso_d = Pers_V2{l}.pers(folds_d{k},:)';
mFD_d = mFD(folds_d{k},:)';


%% for robust regression cross validation
for y = 1:size(perso_d,2);
    left_perso = perso_d(y); %Left out subject personality
    left_disc = discovery(y,:); %Left out subject connectome
    perso_d_cv = perso_d;
    perso_d_cv(y) = []; %remove the left out patient
%     mFD_d_cv = mFD_d;
%     mFD_d_cv(y) = []; %remove motion of the left out patient
    discovery_cv = discovery;
    discovery_cv(:,y) = []; %remove the left out patient

    
for i = 1:size(discovery,1);
    [b,stats] = robustfit(perso_d_cv, discovery_cv(i,:)); %robust regression on all links
%     [b2,stats2] = robustfit(mFD_d_cv, discovery_cv(i,:)); %identify effets of motion
    all_p(i,:) = stats.p;
    all_t(i,:) = stats.t;
%     all_motion_p(i,:) = stats2.p;
%     all_motion_t(i,:) = stats2.t;
end

stat_p = all_p(:,2);
stat_t = all_t(:,2);
% motion_p = all_motion_p(:,2);
% motion_t = all_motion_t(:,2);
idx_connections = find(stat_p < 0.05); %find significant links; to comment for cross validation within sign links
t_sign{y} = stat_t(idx_connections);
p_sign{y} = stat_p(idx_connections);
idx_connec{y} = idx_connections;
% motion_connec{y} = find(motion_p(idx_connections)<0.05);
% idx_connec_fd{y} = idx_connections;
% idx_connec_fd{y}(motion_connec{y})= [];

idx_pos_connections = find(stat_t(idx_connec{y}) > 0); %find positive links
idx_neg_connections = find(stat_t(idx_connec{y}) < 0); %find negative links
idx_connec_pos{y} = idx_pos_connections;
idx_connec_neg{y} = idx_neg_connections;

end

%% identify robust connections across LOOCV

LOOCV_discovery_pos = zeros(size(perso_d,2),size(discovery,1));
LOOCV_discovery_neg = zeros(size(perso_d,2),size(discovery,1));
for g = 1:size(perso_d,2);
    LOOCV_discovery_pos(g,idx_connec{g}(idx_connec_pos{g})) = 1; %% change pos or neg connections
    LOOCV_discovery_neg(g,idx_connec{g}(idx_connec_neg{g})) = 1; %% change pos or neg connections
end

idx_all_pos = sum(LOOCV_discovery_pos,1);
idx_all_neg = sum(LOOCV_discovery_neg,1);


%% identify common connections across LOOCV
idx_pers_pos = mean(LOOCV_discovery_pos,1);
idx_total_pos = find(idx_pers_pos==1);
idx_pers_neg = mean(LOOCV_discovery_neg,1);
idx_total_neg = find(idx_pers_neg==1);

idx_con_pos{l}{k}=idx_total_pos;
idx_con_neg{l}{k}=idx_total_neg;


end
end

