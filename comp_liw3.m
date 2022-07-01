function [ss_liw ps_liw tp_liw pt_liw] = comp_liw3(ss,tp,pp)
if isempty(ss)
    ss = NaN(size(tp));
end

dp = nanmean(nanmean(diff(pp)));
ss_tmp = ss;tp_tmp = tp;gg = sw_dens(ss,tp,zeros(size(tp)))-1000;
ss_tmp(gg<28.95 | gg>29.115) = NaN;ss_tmp = median_av2(ss_tmp,round(5/dp)); %%%% en chaque z d'un profil, on calcul la médiane des index valeurs autour de la valeur à ce z
tp_tmp(gg<28.95 | gg>29.115) = NaN;tp_tmp = median_av2(tp_tmp,round(5/dp));


ss_tmp(:,find(nanmax(ss_tmp)>38.9)) = NaN;tp_tmp(:,find(nanmax(ss_tmp)>38.9)) = NaN;
pp_tmp = pp; pp_tmp(isnan(tp)) = NaN;
%pp_tmp = pp; pp_tmp(unique(union(isnan(ss)) = NaN;
zmax = nanmax(pp_tmp);
% iqnan = find(abs(diff([NaN(1,length(ss(1,:))) ; ss_tmp]))>0.007);
% ss_tmp(iqnan) = NaN;tp_tmp(iqnan) = NaN;
iq = find(zmax<350 | nanmax(gg)<28.95);
ss_tmp(:,iq) = NaN;tp_tmp(:,iq) = NaN;


%[tp_liw pt_liw] = nanmax(tp_tmp(round(100/dp):end,:));
[ss_liw ps_liw] = nanmax(ss_tmp(round(100/dp):end,:));
%tp_liw(pt_liw==1)=NaN;pt_liw(pt_liw==1)=NaN;pt_liw=round(100/dp)+pt_liw*dp;
ss_liw(ps_liw==1)=NaN;ps_liw(ps_liw==1)=NaN;ps_liw=round(100/dp)+ps_liw*dp;
tp_tmp(1:100,:)=NaN;
tp_tmp(1000:end,:)=NaN;
 
%%%movng average of the temperature matrice
for d=1:size(tp_tmp,2)
tp_tmp2(:,d) = mean_av(tp_tmp(:,d),round(5/dp)); %%%% en chaque z d'un profil, on calcul la moyenne des index valeurs autour de la valeur à ce z
end

%tp_tmp2=tp_tmp;
%%% use the smoothed matrice to calculate the partial slopes of each
%%% profile
for pro=1:size(tp_tmp2,2)
    vert=1;
    for start_ind=100:50:550
        if size(tp_tmp2,1)>start_ind+100
        last_ind=start_ind+100;
        else
            last_ind=size(tp_tmp2,1);
        end
        %faire une regression linéaire des points de tp_tmp entre la profondeur start_ind et last_ind
        indi= [start_ind:last_ind];
        ind= indi(find(~isnan(tp_tmp2(indi,pro))));
        P=polyfit(pp_tmp(ind,pro),tp_tmp2(ind,pro),1);
        warning('off','all')
        slope(vert,pro)=P(1);
        vert=vert+1;
    end
end


for pro=1:size(tp_tmp2,2)
    [posi,~]=find(slope(:,pro)>0,1,'last'); %%find the last 
    [nega,~]=find(slope(:,pro)<0,1,'last');
    if ~isempty(posi)
        ind_sl_pos(pro)=posi;
    else
        ind_sl_pos(pro)=NaN;
    end
    if ~isempty(nega)
        ind_sl_neg(pro)=nega;
    else
        ind_sl_neg(pro)=NaN;
    end
end

to_keep=find(ind_sl_pos-ind_sl_neg<0);
to_erase=setdiff([1:size(tp_tmp,2)],to_keep);
id_sl_pos(to_erase)=NaN;
tp_tmp(:,to_erase)=NaN;
tt=[100:50:550];

for g=1:size(tp_tmp,2)
    if ~isnan(ind_sl_pos(g))
    tp_tmp(1:tt(ind_sl_pos(g)),g)=NaN;
    else
        continue
    end
end
[tp_liw pt_liw] = nanmax(tp_tmp);
tp_liw(pt_liw==1)=NaN;pt_liw(pt_liw==1)=NaN;



%%% if sal_liw detected but no tpot_liw look for slopes in smaller trasects



intersect_ind= setdiff([find(isnan(pt_liw))],[intersect(find(isnan(pt_liw)),find(isnan(ss_liw)))]);
tp_tmp3=tp_tmp2(:,intersect_ind);
for i=1:length(intersect_ind)
tp_tmp3(ps_liw(intersect_ind(i)):end)=NaN;
end
if ~isempty (intersect_ind)
    
    
for pro=1:size(tp_tmp3,2)
    vert=1;
    for start_ind=100:10:540
        if size(tp_tmp2,1)>start_ind+30
        last_ind=start_ind+30;
        else
            last_ind=size(tp_tmp3,1);
        end
        %faire une regression linéaire des points de tp_tmp entre la profondeur start_ind et last_ind
        indi= [start_ind:last_ind];
        ind= indi(find(~isnan(tp_tmp3(indi,pro))));
        P=polyfit(pp_tmp(ind,pro),tp_tmp3(ind,pro),1);
        warning('off','all')
        slope(vert,pro)=P(1);
        vert=vert+1;
    end
end
ind_sl_pos=[];
ind_sl_neg=[];

for pro=1:size(tp_tmp3,2)
    [posi,~]=find(slope(:,pro)>0,1,'last'); %%find the last 
    [nega,~]=find(slope(:,pro)<0,1,'last');
    if ~isempty(posi)
        ind_sl_pos(pro)=posi;
    else
        ind_sl_pos(pro)=NaN;
    end
    if ~isempty(nega)
        ind_sl_neg(pro)=nega;
    else
        ind_sl_neg(pro)=NaN;
    end
end

to_keep=find(ind_sl_pos-ind_sl_neg<0);
to_erase=setdiff([1:size(tp_tmp3,2)],to_keep);
ind_sl_pos(to_erase)=NaN;
tp_tmp3(:,to_erase)=NaN;
tt=[100:10:540];

for g=1:size(tp_tmp3,2)
    if ~isnan(ind_sl_pos(g))
    tp_tmp3(1:tt(ind_sl_pos(g)),g)=NaN;
    else
        continue
    end
end
[tp_liw3 pt_liw3] = nanmax(tp_tmp3);
tp_liw3(pt_liw3==1)=NaN;pt_liw3(pt_liw3==1)=NaN;
    
tp_liw(intersect_ind)=tp_liw3;
pt_liw(intersect_ind)=pt_liw3;

end    
 







