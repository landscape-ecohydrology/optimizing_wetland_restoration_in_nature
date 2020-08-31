clear
clc

SIM='10';

filename=['WETLAND_',SIM,'_pareto.mat'];                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
load (filename)

ID=round(LAND_data(:,1)/1000);
CROPLAND_area_random=LAND_data(:,2);
CROPLAND_area_no_ag=LAND_data(:,5);
CROPLAND_area_targeted=LAND_data(:,8);
PASTLAND_area_random=LAND_data(:,3);
PASTLAND_area_no_ag=LAND_data(:,6);
PASTLAND_area_targeted=LAND_data(:,9);
OTHER_area_random=LAND_data(:,4);
OTHER_area_no_ag=LAND_data(:, 7);
OTHER_area_targeted=LAND_data(:,10);

% RANDOM_area=CROPLAND_area_random+PASTLAND_area_random+OTHER_area_random;
% CO_LO_area=CROPLAND_area_targeted+PASTLAND_area_targeted+OTHER_area_targeted;

% RANDOM_area_total=sum(RANDOM_area);
% CO_LO_area_total=sum(CO_LO_area);

A=[CROPLAND_area_random CROPLAND_area_no_ag CROPLAND_area_targeted];
B=[PASTLAND_area_random PASTLAND_area_no_ag PASTLAND_area_targeted];
C=[OTHER_area_random OTHER_area_no_ag OTHER_area_targeted];
WETLAND_area=cat(3,A,B,C);
                                       
filename='cropland_rents_state_2017.csv';
T=readtable(filename);
STATE_ID=T.StateANSI;
filename='pastland_rents_state_2017.csv';
CROP_rent=T.Value;
T2=readtable(filename);
STATE_ID2=T2.StateANSI;
PAST_rent=T2.Value;

% Cropland Rents

% STATE_ID=NaN(size(T,1),1);
% CROP_rent=NaN(size(T,1),1);
% 
% for i=1:size(T,1)
%     STATE_ID(i,1)=str2double(string(cell2mat(ID_temp(i))));
%     CROP_rent(i,1)=str2double(string(cell2mat(CROP_rent_temp(i)))) / 0.404686;  % $ per hectare
% end

% Pastureland Rents

% STATE_ID2=NaN(size(T2,1),1);
% PAST_rent=NaN(size(T2,1),1);
% 
% for i=1:size(T2,1)
%     STATE_ID2(i,1)=str2double(string(cell2mat(ID_temp2(i))));
%     PAST_rent(i,1)=str2double(string(cell2mat(PAST_rent_temp(i)))) / 0.404686;  % $ per hectare;
% end

PUBLIC_rent=zeros(size(T2,1),1);

RENT_cost=zeros(size(LAND_data,1),3,3);
RENT_buffer_cost=zeros(size(LAND_data,1),3,3);
DESIGN_cost=zeros(size(LAND_data,1),2,2);
CONSTRUCT_cost=zeros(size(LAND_data,1),2,2);
PLANT_cost=zeros(size(LAND_data,1),2,2);
SEED_cost=zeros(size(LAND_data,1),2,2);
SEEDING_cost=zeros(size(LAND_data,1),2,2);
WEIR_cost=zeros(size(LAND_data,1),2,2);
CONTROL_cost=zeros(size(LAND_data,1),2,2);
TIME_cost=zeros(size(LAND_data,1),2,2);
REPLACE_gate_cost=zeros(size(LAND_data,1),2,2);
REPLACE_control_cost=zeros(size(LAND_data,1),2,2);
RENT_all=zeros(size(LAND_data,1),5,3);

BUFFER=0.035; % percent of wetland area
DESIGN=1000 / 0.404686 / 50; % USD per wetland hectare per year
CONSTRUCT= 1500 / 0.404686 / 50; % USD per wetland hectare per year
PLANT= 640  / 0.404686 / 50; % USD per wetland hectare per year
SEED = 131  / 0.404686 / 50; % USD per hectare per year
SEEDING = 40 / 0.404686 / 50; % USD per hectare per year
WEIR = 600 / 0.404686 / 50; % USD per hectare per year
CONTROL = 2100 / 0.404686 / 50; % USD per hectare per year
TIME = 3.06 / 0.404686 / 47; % USD per hectare per year
REPLACE_gate = 15 / 0.404686 * 6 / 50; % USD per hectare per year
REPLACE_control = 935 / 0.404686 / 50; % USD per hectare per year

COST_sum=zeros(size(LAND_data,1),2,2);

ent=size(LAND_data,1);

for i=1:ent
	for j=1:3 % scenario
        for k=1:3 % landuse
            if k==1
                idx=STATE_ID==ID(i);
                RENT=CROP_rent(idx);
            elseif k==2
                idx=STATE_ID2==ID(i);
                RENT=PAST_rent(idx);
            elseif k==3
                RENT=0;
            end
            AREA=WETLAND_area(i,j,k);
            if AREA>0
                RENT_cost(i,j,k)=AREA*RENT;  
                RENT_buffer_cost(i,j,k)=AREA*BUFFER*RENT;
                RENT_all(i,j,k)=RENT_cost(i,j,k)+RENT_buffer_cost(i,j,k);
                DESIGN_cost(i,j,k)=AREA*DESIGN;  
                CONSTRUCT_cost(i,j,k)=AREA*CONSTRUCT;
                PLANT_cost(i,j,k)=AREA*PLANT;  
                SEED_cost(i,j,k)=AREA*SEED;
                SEEDING_cost(i,j,k)=AREA*SEEDING;
                WEIR_cost(i,j,k)=AREA*WEIR;
                CONTROL_cost(i,j,k)=AREA*CONTROL;
                TIME_cost(i,j,k)=AREA*TIME;  
                REPLACE_gate_cost(i,j,k)=AREA*REPLACE_gate;
                REPLACE_control_cost(i,j,k)=AREA*REPLACE_control;
            end
        end
    end
end

RENT_all_total=NaN(3,1);
DESIGN_total=NaN(3,1);
CONSTRUCT_total=NaN(3,1);
PLANT_total=NaN(3,1);
SEED_total=NaN(3,1);
SEEDING_total=NaN(3,1);
WEIR_total=NaN(3,1);
CONTROL_total=NaN(3,1);
TIME_total=NaN(3,1);
REPLACE_gate_total=NaN(3,1);
REPLACE_control_total=NaN(3,1);
ALL_total=NaN(3,1);
per_RENT_total=NaN(3,1);

for j=1:3  % 3 scenarios
	RENT_all_total(j)=sum(RENT_all(:,j,:),'all');
    DESIGN_total(j)=sum(DESIGN_cost(:,j,:),'all');
    CONSTRUCT_total(j)=sum(CONSTRUCT_cost(:,j,:),'all');
    PLANT_total(j)=sum(PLANT_cost(:,j,:),'all');
    SEED_total(j)=sum(SEED_cost(:,j,:),'all');
    SEEDING_total(j)=sum(SEEDING_cost(:,j,:),'all');
    WEIR_total(j)=sum(WEIR_cost(:,j,:),'all');
    CONTROL_total(j)=sum(CONTROL_cost(:,j,:),'all');
    TIME_total(j)=sum(TIME_cost(:,j,:),'all');
    REPLACE_gate_total(j)=sum(REPLACE_gate_cost(:,j,:),'all');
    REPLACE_control_total(j)=sum(REPLACE_control_cost(:,j,:),'all');
    ALL_total(j)=RENT_all_total(j)+DESIGN_total(j)+CONSTRUCT_total(j)+PLANT_total(j)+SEED_total(j)+SEEDING_total(j)+WEIR_total(j)+CONTROL_total(j)+TIME_total(j)+REPLACE_gate_total(j)+REPLACE_control_total(j);
    per_RENT_total(j)=RENT_all_total(j)/ALL_total(j);
end

save (['PARETO_results/PARETO_COST_results_',SIM,'.mat'])

% CROPLAND_cost_total_random=sum(CROPLAND_cost_random);
% CROPLAND_cost_total_targeted=sum(CROPLAND_cost_targeted);
% CROPLAND_lost_total_targeted=nansum(CROPLAND_area_targeted);
% 
% 
% PASTLAND_cost_total_random=sum(PASTLAND_cost_random);
% PASTLAND_cost_total_targeted=sum(PASTLAND_cost_targeted);
% PASTLAND_lost_total_targeted=nansum(PASTLAND_area_targeted);
% 
% 
% AG_cost_total_random=PASTLAND_cost_total_random+CROPLAND_cost_total_random;
% AG_cost_total_targeted=PASTLAND_cost_total_targeted+CROPLAND_cost_total_targeted;
% AG_lost_total_targeted=PASTLAND_lost_total_targeted+CROPLAND_lost_total_targeted;
