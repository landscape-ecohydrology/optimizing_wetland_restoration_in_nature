clear
clc
close all

NAME='PARETO_FIG5_RESULTS/0622_005'
TOT=5.1E11; % total wetland area
INCR=0.005;
max_syst_wetland_area=TOT*INCR; % 10% increase in wetland area 

tic

% initiate monte carlo 

ITER=5; % number of iterations
np=6; % number of parameters to vary in monte carlo

LH=lhsdesign(ITER,np);

% What fraction of N surplus reaches the wetland?

NS_frac_o=0.4;

% How many watersheds?

n_wshd=3262;  % number of watersheds
p_wshd=97;   % number of pixels per watershed
p_area=25; %km2
p_wshd_vec=1:p_wshd;

% Watershed Size

size_wshd=p_wshd*p_area; % km2

% Import Nitrogen Surplus Data

load ('US_data_v3_no_shrub.mat','M2')

M2_append=zeros(3,size(M2,2));
M2=[M2;M2_append];

NS=reshape(M2(:,3),p_wshd,n_wshd);
NS_original=NS;
LU=reshape(M2(:,4),p_wshd,n_wshd);
LU_open=reshape(M2(:,5),p_wshd,n_wshd);

% Generate wetland size distribution and related parameters

mu=500;    % mean wetland size
b=1.84;     % power law exponent
a=mu*(b-1)/b;  % power law coefficient
n=10^6; % number of wetlands generated

SA = randraw('pareto', [a, b], [1 n]); % wetland surface areas
idx=SA<=10^5;
SA=SA(idx);

a_SA_o=1.51;
b_SA_o=0.23;

a_tau_o=0.38;
b_tau_o=-0.91;

CA_SA_ratio_o=8.4; % range 3-20

load '/Volumes/GoogleDrive/My Drive/Research/Matlab/US Wetlands/Catchment Area/CA_trunc.mat'; % pull in values for CA distribution

% Create Monte Carlo Parameter Set

NS_frac_unc=0.1;
a_SA_unc=0.07;
b_SA_unc=0.02;
a_tau_unc=0.07;
b_tau_unc=0.05;
CA_ratio_LB=3;
CA_ratio_UB=20;


for mc=1:ITER

    NS_frac(mc,1)=(NS_frac_o-NS_frac_unc)+(2*(NS_frac_unc)*LH(mc,1));
    a_SA(mc,1)=(a_SA_o-a_SA_unc)+(2*(a_SA_unc)*LH(mc,2));
    b_SA(mc,1)=(b_SA_o-b_SA_unc)+(2*(b_SA_unc)*LH(mc,3));
    a_tau(mc,1)=(a_tau_o-a_tau_unc)+(2*(a_tau_unc)*LH(mc,4));
    b_tau(mc,1)=(b_tau_o-b_tau_unc)+(2*(b_tau_unc)*LH(mc,5));
    CA_SA_ratio(mc,1)=WET_trunc(round(length(WET_trunc)*LH(mc,6)));
    %CA_SA_ratio(mc,1)=8.4;
    %CA_SA_ratio(mc,1)=CA_ratio_LB+(CA_ratio_UB-CA_ratio_LB)*LH(mc,6);

tau(mc,:) = a_SA(mc)*(SA.^b_SA(mc)); % wetland residence time, days
k_rate(mc,:)=a_tau(mc)*(tau(mc,:).^(b_tau(mc))); % removal rate constant
%CA = 10.^(0.8 * log10(SA) + 2); % equation from Fred/Lake CAT 
CA(mc,:)=CA_SA_ratio(mc)*SA;  % catchment area, use geometric mean from Wu & Lane

end


%%%%%%%%%%%%%%%%%%%%%%%%
% Place New Wetlands %%%
%%%%%%%%%%%%%%%%%%%%%%%%

n_sim=3

% Define Variables

N_rem_total_wshd=NaN(n_wshd,ITER,n_sim); % tons N
N_surp_total_wshd=NaN(n_wshd,n_sim);  % tons N
N_per_removal_wshd=NaN(n_wshd,ITER,n_sim); % fraction
C_area_total_wshd=NaN(n_wshd,ITER,n_sim);
C_area_fraction=NaN(n_wshd,ITER,n_sim);
W_num_total_wshd=NaN(n_wshd,ITER,n_sim);

N_rem_total=NaN(n_sim,ITER);
%W_area_total=NaN(n_sim,ITER);
%W_num_total=NaN(n_sim,ITER);

CROPLAND_lost_random_mc=zeros(p_wshd,n_wshd,ITER);
PASTLAND_lost_random_mc=zeros(p_wshd,n_wshd,ITER);
OTHER_lost_random_mc=zeros(p_wshd,n_wshd,ITER);
        
CROPLAND_lost_co_location_mc=zeros(p_wshd,n_wshd,ITER);
PASTLAND_lost_co_location_mc=zeros(p_wshd,n_wshd,ITER);
OTHER_lost_co_location_mc=zeros(p_wshd,n_wshd,ITER);

CROPLAND_lost=zeros(p_wshd,n_wshd,n_sim);
PASTLAND_lost=zeros(p_wshd,n_wshd,n_sim);
OTHER_lost=zeros(p_wshd,n_wshd,n_sim);


for mc=1:ITER
    
    CROPLAND_lost=zeros(p_wshd,n_wshd,n_sim);
    PASTLAND_lost=zeros(p_wshd,n_wshd,n_sim);
    OTHER_lost=zeros(p_wshd,n_wshd,n_sim);
    
    CROPLAND_C_area=zeros(p_wshd,n_wshd,n_sim,ITER);
    PASTLAND_C_area=zeros(p_wshd,n_wshd,n_sim,ITER);
    OTHER_C_area=zeros(p_wshd,n_wshd,n_sim,ITER);


    for i=1:3
        
        i
    
    % 1 = random, 2 = no cropland loss, 3 = co-location
    
    NS=NS_original;
    
    if i<=2
        nest_factor=1;
    else
        nest_factor=1;
    end
    
    %max_syst_wetland_area=1.7978E11;  
    %max_syst_wetland_area=1.115E11; % 25% increase in wetland area  
    %max_syst_wetland_area=5.1E11; % Fred's current NWI wetland area  
    %max_syst_wetland_area=1.3E9; % CRP Remaining Wetland Area #0101a
    %max_syst_wetland_area=2.5E10; % 5% increase in wetland area #0101b
    %max_syst_wetland_area=5.1E9; % 1% increase in wetland area #0101b
    %max_syst_wetland_area=1.02E10; % 2% increase in wetland area 
    %max_syst_wetland_area=2.55E9; % 0.5% increase in wetland area #0102a
    %max_syst_wetland_area=3.825E10; % 7.5% increase in wetland area 
    %max_syst_wetland_area=5.1E10; % 10% increase in wetland area  
    
    
    max_wetland_catch_pixel_area=nest_factor*p_area*1000^2; %m2

    W_area_wshd=zeros(p_wshd,n_wshd);
    N_rem_wshd=zeros(p_wshd,n_wshd);
    C_area_wshd=zeros(p_wshd,n_wshd);
    W_num_wshd=zeros(p_wshd,n_wshd);
    
    
    NS_test=reshape(NS,numel(NS),1);
    idx=NS_test>0;
    NS_test=NS_test(idx);
    
    if i<=2
        NS_thresh=0;
    else
        thresh_val=90;
        NS_thresh=prctile(NS_test,thresh_val);
        thresh_progress=0.03*mu*n_wshd;
    end
    
    
    W_progress=0;

    while sum(W_area_wshd(:,:),'all')<0.999*max_syst_wetland_area
        for j=1:n_wshd
            wshd_pixel=randperm(numel(p_wshd_vec), 1);
             if LU_open(wshd_pixel,j)==1
                if NS(wshd_pixel,j)>=NS_thresh
                    GO=1;
                    if i==2
                        if LU(wshd_pixel,j)>0
                            GO=0;
                        end
                    end
                    if GO==1
                        rnum=round((length(SA)-1).*rand+1);
                        W_area=SA(rnum);
                        W_tau=tau(mc,rnum);
                        W_k=k_rate(mc,rnum);
                        C_area=CA(mc,rnum);
                        C_area_pixel_temp=C_area_wshd(wshd_pixel,j)+C_area; %Catchment Area
                        if C_area_pixel_temp<max_wetland_catch_pixel_area
                            W_num_wshd(wshd_pixel,j)=W_num_wshd(wshd_pixel,j)+1;
                            W_area_wshd(wshd_pixel,j)=W_area_wshd(wshd_pixel,j)+W_area;;
                            C_area_wshd(wshd_pixel,j)=C_area_wshd(wshd_pixel,j)+C_area;
                            N_surp=NS_frac(mc)*C_area/10^4*NS(wshd_pixel,j); %kg                   
                            N_rem_wshd(wshd_pixel,j)=N_surp*(1-exp(-W_k*W_tau))+N_rem_wshd(wshd_pixel,j); %kg
                            %NS(wshd_pixel,j)=((NS(wshd_pixel,j)*p_area*1000^2/10^4)-N_surp)/(p_area*1000^2/10^4); %remaining NS, kg/ha
                            if LU(wshd_pixel,j)==1
                                PASTLAND_lost(wshd_pixel,j,i)=PASTLAND_lost(wshd_pixel,j,i)+W_area;
                                PASTLAND_C_area(wshd_pixel,j,i)=PASTLAND_C_area(wshd_pixel,j,i)+C_area;
                            elseif LU(wshd_pixel,j)==2
                                CROPLAND_lost(wshd_pixel,j,i)=CROPLAND_lost(wshd_pixel,j,i)+W_area;
                                CROPLAND_C_area(wshd_pixel,j,i)=CROPLAND_C_area(wshd_pixel,j,i)+C_area;
                            elseif LU(wshd_pixel,j)==0
                                OTHER_lost(wshd_pixel,j,i)=OTHER_lost(wshd_pixel,j,i)+W_area;
                                OTHER_C_area(wshd_pixel,j,i)=OTHER_C_area(wshd_pixel,j,i)+C_area;
                            end
                        end
                    end
                end
            end       
        end
        if i==3
            if (sum(W_area_wshd,'all')-W_progress)<thresh_progress
                thresh_val=thresh_val-1;
                NS_thresh=prctile(NS_test,thresh_val);
                if thresh_val<=0
                    thresh_val=1;
                    NS_thresh=prctile(NS_test,thresh_val);
                end
            end
            W_progress=sum(W_area_wshd,'all');  
        end
    end
    
    
	for j=1:n_wshd
        N_rem_total_wshd(j,mc,i)=sum(N_rem_wshd(:,j))/1000;         % tons N
    end
% % %             N_surp_total_wshd(j,i)=sum(NS_original(:,j)*size_wshd*1000^2/10^4/1000); % tons N
% % %             N_per_removal_wshd(j,mc,i)=N_rem_total_wshd(j,i)/N_surp_total_wshd(j,i); % fraction
% % %             C_area_total_wshd(j,mc,i)=sum(C_area_wshd(:,j));
% % %             C_area_fraction(j,mc,i)=C_area_total_wshd(j,mc,i)/(size_wshd*1000^2);
% % %             W_num_total_wshd(j,mc,i)=sum(W_num_wshd(:,j));
% % %     end
    
    if i==1
        N_rem_random_mc(:,:,mc)=N_rem_wshd(:,:);
        C_area_wshd_random_mc(:,:,mc)=C_area_wshd(:,:);
        W_area_wshd_random_mc(:,:,mc)=W_area_wshd(:,:);
        W_num_wshd_random_mc(:,:,mc)=W_num_wshd(:,:);
        
        CROPLAND_lost_random_mc(:,:,mc)=CROPLAND_lost(:,:,i);
        PASTLAND_lost_random_mc(:,:,mc)=PASTLAND_lost(:,:,i);
        OTHER_lost_random_mc(:,:,mc)=OTHER_lost(:,:,i);
        
        CROPLAND_C_area_random_mc(:,:,mc)=PASTLAND_C_area(:,:,i);
        PASTLAND_C_area_random_mc(:,:,mc)=PASTLAND_C_area(:,:,i);
        PASTLAND_C_area_random_mc(:,:,mc)=PASTLAND_C_area(:,:,i);

    
    elseif i==2
        
        N_rem_no_crop_mc(:,:,mc)=N_rem_wshd(:,:);
        C_area_wshd_no_crop_mc(:,:,mc)=C_area_wshd(:,:);
        W_area_wshd_no_crop_mc(:,:,mc)=W_area_wshd(:,:);
        W_num_wshd_no_crop_mc(:,:,mc)=W_num_wshd(:,:);

    elseif i==3
        
         N_rem_co_location_mc(:,:,mc)=N_rem_wshd(:,:); 
        C_area_wshd_co_location_mc(:,:,mc)=C_area_wshd(:,:);
        W_area_wshd_co_location_mc(:,:,mc)=W_area_wshd(:,:);
        W_num_wshd_co_location_mc(:,:,mc)=W_num_wshd(:,:);
        
        CROPLAND_lost_co_location_mc(:,:,mc)=CROPLAND_lost(:,:,i);
        PASTLAND_lost_co_location_mc(:,:,mc)=PASTLAND_lost(:,:,i);
        OTHER_lost_co_location_mc(:,:,mc)=OTHER_lost(:,:,i);
        
        
        CROPLAND_C_area_co_location_mc(:,:,mc)=PASTLAND_C_area(:,:,i);
        PASTLAND_C_area_co_location_mc(:,:,mc)=PASTLAND_C_area(:,:,i);
        PASTLAND_C_area_co_location_mc(:,:,mc)=PASTLAND_C_area(:,:,i);
        
    end

	N_rem_total(i,mc)=sum(squeeze(N_rem_total_wshd(:,mc,i))); % tons


    end


toc


% % % 
% % % N_surp_total=sum(NS_original,'all')*p_area*10^6/10^4; % kg
% % % N_surp_total_tons=N_surp_total/1000; % tons
% % % N_per_removal(i,mc)=N_rem_total(i,mc)./N_surp_total_tons; %fraction



% % % figure(8)
% % % bar([N_rem_total(2) N_rem_total(1) N_rem_total(3)])
% % % ylim([0 10^6])

% % % filename=['FIGURES/N_rem_bar_',NAME,'.jpg'];
% % % saveas(figure(8),filename);
% % % 


% % % figure(9)
% % % bar([N_per_removal(2) N_per_removal(1) N_per_removal(3)])
% % % 
% % % filename=['FIGURES/N_per_removal_bar_',NAME,'.jpg'];
% % % saveas(figure(9),filename);
% % % 
% % % 

% % % for i=1:n_sim
% % %     CROPLAND_lost_total(mc,i)=sum(squeeze(CROPLAND_lost(:,:,i)),'all')/10^4; %ha
% % %     PASTLAND_lost_total(mc,i)=sum(squeeze(PASTLAND_lost(:,:,i)),'all')/10^4; %ha
% % %     OTHER_lost_total(mc,i)=sum(squeeze(OTHER_lost(:,:,i)),'all')/10^4; %ha
% % % end

end

N_rem_mean=mean(N_rem_total,2)/10^3;  % in %Tg
N_rem_std=std(N_rem_total,0,2)/10^3;  % in Tg

% % % figure(10)
% % % subplot(2,1,1)
% % % bar([CROPLAND_lost_total(1) CROPLAND_lost_total(3)])
% % % subplot(2,1,2)
% % % bar([PASTLAND_lost_total(1) PASTLAND_lost_total(3)])
% % % 
% % % filename=['FIGURES/AGLAND_lost_',NAME,'.jpg'];
% % % saveas(figure(10),filename);

% Get mean values

%N_rem_random=nanmean(N_rem_random_mc,3);
%N_rem_no_crop=nanmean(N_rem_no_crop_mc,3);
N_rem_co_location=nanmean(N_rem_co_location_mc,3);

save ([NAME,'MONTE_CARLO_3_runs.mat'])

%N_rem_random_std=nanstd(N_rem_random_mc,3);
%N_rem_no_crop_std=nanstd(N_rem_no_crop_mc,3);
N_rem_co_location_std=nanstd(N_rem_co_location_mc,0,3);

%W_area_wshd_random=nanmean(W_area_wshd_random_mc,3);
%W_area_wshd_no_crop=nanmean(W_area_no_crop_mc,3);
W_area_wshd_co_location=nanmean(W_area_co_location,3);

%W_area_wshd_random_std=nanstd(W_area_wshd_random_mc,3);
%W_area_wshd_no_crop_std=nanstd(W_area_no_crop_mc,3);
W_area_wshd_co_location_std=nanstd(W_area_co_location,3);

%CROPLAND_lost_random=nanmean(CROPLAND_lost_random_mc,3);
%PASTLAND_lost_random=nanmean(PASTLAND_lost_random_mc,3);
OTHER_lost_random=nanmean(OTHER_lost_random_mc,3);

%CROPLAND_lost_random_std=nanstd(CROPLAND_lost_random_mc,3);
%PASTLAND_lost_random_std=nanstd(PASTLAND_lost_random_mc,3);
OTHER_lost_random_std=nanstd(OTHER_lost_random_mc,3);

%CROPLAND_C_area_co_location=nanmean(CROPLAND_C_area_co_location_mc,3);
%PASTLAND_C_area_co_location=nanmean(PASTLAND_C_area_co_location_mc,3);
OTHER_C_area_co_location=nanmean(OTHER_C_area_co_location_mc,3);

%CROPLAND_C_area_co_location_std=nanstd(CROPLAND_C_area_co_location_mc,3);
%PASTLAND_C_area_co_location_std=nanstd(PASTLAND_C_area_co_location_mc,3);
OTHER_C_area_co_location_std=nanstd(OTHER_C_area_co_location_mc,3);

%N_rem_random=reshape(N_rem_random,size(M2,1),1);
%N_rem_no_crop=reshape(N_rem_no_crop,size(M2,1),1);
N_rem_co_location=reshape(N_rem_co_location,size(M2,1),1);

%N_rem_random_std=reshape(N_rem_random_std,size(M2,1),1);
%N_rem_no_crop_std=reshape(N_rem_no_crop_std,size(M2,1),1);
N_rem_co_location_std=reshape(N_rem_co_location_std,size(M2,1),1);

%W_area_random=reshape(W_area_wshd_random,size(M2,1),1);
%W_area_no_crop=reshape(W_area_wshd_no_crop,size(M2,1),1);
W_area_co_location=reshape(W_area_wshd_co_location,size(M2,1),1);

%W_area_random_std=reshape(W_area_random_std,size(M2,1),1);
%W_area_no_crop_std=reshape(W_area_no_crop_std,size(M2,1),1);
W_area_co_location_std=reshape(W_area_co_location_std,size(M2,1),1);

% % % CROPLAND_lost_random=reshape(CROPLAND_lost(:,:,1),size(M2,1),1)/10^4; % ha
% % % CROPLAND_lost_co_location=reshape(CROPLAND_lost(:,:,3),size(M2,1),1)/10^4; % ha
% % % CROPLAND_C_area_random=reshape(CROPLAND_C_area(:,:,1),size(M2,1),1)/10^4; % ha
% % % CROPLAND_C_area_co_location=reshape(CROPLAND_C_area(:,:,3),size(M2,1),1)/10^4; % ha
% % % PASTLAND_lost_random=reshape(PASTLAND_lost(:,:,1),size(M2,1),1)/10^4; % ha
% % % PASTLAND_lost_co_location=reshape(PASTLAND_lost(:,:,3),size(M2,1),1)/10^4; % ha
% % % PASTLAND_C_area_random=reshape(PASTLAND_C_area(:,:,1),size(M2,1),1)/10^4; % ha
% % % PASTLAND_C_area_co_location=reshape(PASTLAND_C_area(:,:,3),size(M2,1),1)/10^4; % ha
% % % OTHER_lost_random=reshape(OTHER_lost(:,:,1),size(M2,1),1)/10^4; % ha
% % % OTHER_lost_co_location=reshape(OTHER_lost(:,:,3),size(M2,1),1)/10^4; % ha
% % % OTHER_C_area_random=reshape(OTHER_C_area(:,:,1),size(M2,1),1)/10^4; % ha
% % % OTHER_C_area_co_location=reshape(OTHER_C_area(:,:,3),size(M2,1),1)/10^4; % ha

%W_per_random=W_area_random/(25*1000^2);
%W_per_no_crop=W_area_no_crop/(25*1000^2);
W_per_co_location=W_area_co_location/(25*1000^2);

M2=[M2 N_rem_co_location W_area_co_location W_per_co_location];
AGLAND_lost_data=[M2 CROPLAND_lost_co_location PASTLAND_lost_co_location OTHER_lost_co_location ];

STATE_geoid=floor(M2(:,2)/10^3);

M2=[M2 STATE_geoid];

filename=['M2_',NAME,'.csv'];
writetable(table(M2),filename)

save (['CROPLAND_data_mc',NAME,'.mat'],'AGLAND_lost_data')
save (['/Volumes/GoogleDrive/My Drive/Research/Matlab/US Wetlands/Cropland Cost Calculations/CROPLAND_data_',NAME,'.mat'],'AGLAND_lost_data')

save (['WETLANDS_US_v4_mc',NAME,'.mat'])
save (['/Volumes/GoogleDrive/My Drive/Research/Matlab/US Wetlands/landuse calculations/WETLANDS_US_v4_',NAME,'.mat'])




