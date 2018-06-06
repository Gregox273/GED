h = xtreme_get_rsq_header('../raw_ct_data/raw_data_b.rsq');
mid = (h.fan_scans-2*h.skip_scans+1)/2;
x = [1:h.fan_scans-2*h.skip_scans];
y = abs(mid - x);
%y = sqrt((y.^2 * h.scale^2)/(h.radius^2) + 1);
y = sqrt((h.radius^2)./((y.^2 * h.scale^2) + h.radius^2));
plot(100*(1-[y,y,y]));