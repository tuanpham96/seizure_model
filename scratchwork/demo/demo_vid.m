tmaxviz = 3; %s 

ind_est = data.Rate.time <= tmaxviz; 
est_time = data.Rate.time(ind_est); 
est_rate = data.Rate.rate(:,ind_est); 

ind_stm = data.I_stim.time <= tmaxviz;  
Tsec = data.I_stim.time(ind_stm); 
Istim = data.I_stim.stim(ind_stm); 

stim_ind = data.Stim_Ind;

vid_name = 'demo_vid'; 
% net_arch.create_activity_movie(vid_name, est_time, est_rate, Tsec, Istim, stim_ind); 
%%
x_coord = net_arch.coord.x; 
y_coord = net_arch.coord.y; 

% num_neighs = 10; 
% cnt_nrn = net_arch.centered_neuron.ind; 
% cntpop = net_arch.pick_random_neighbors(num_neighs, cnt_nrn); 
% cntpop = [cntpop, cnt_nrn]; 
% 
% side_nrn = net_arch.side_neuron.ind; 
% sidepop = net_arch.pick_random_neighbors(num_neighs, side_nrn);
% sidepop = [sidepop, side_nrn]; 
% 
% mid_x = (net_arch.centered_neuron.x + net_arch.side_neuron.x)/2;
% mid_y = net_arch.side_neuron.y;
% dist2mid = (x_coord-mid_x).^2 + (y_coord-mid_y).^2; 
% mid_nrn = find(dist2mid == min(dist2mid),1); 
% midpop = net_arch.pick_random_neighbors(num_neighs, mid_nrn); 
% midpop = [midpop, mid_nrn];
% 
% viz_pop = [cntpop, sidepop, midpop];

num_neighs = 20; 
cnt_nrn = net_arch.centered_neuron.ind; 
cntpop = net_arch.pick_random_neighbors(num_neighs, cnt_nrn); 
cntpop = [cntpop, cnt_nrn]; 

unq_x = unique(x_coord); 
cnt_y = net_arch.side_neuron.y;
pop_line = cell(length(unq_x), 1); 
for i = 1:length(unq_x)
    xi = unq_x(i);
    dist2mid = (x_coord-xi).^2 + (y_coord-cnt_y).^2;
    ind_i = find(dist2mid == min(dist2mid),1);
    pop_i = net_arch.pick_random_neighbors(num_neighs, ind_i);
    pop_line{i} = [pop_i, ind_i];
end
viz_pop = [pop_line{:}]; 
%%
x_coord = round(2*(x_coord - min(x_coord))) + 1; 
y_coord = round(2*(y_coord - min(y_coord))) + 1; 

actv_mat = zeros(max(x_coord), max(y_coord), length(est_time));
for i = 1:length(x_coord) 
    actv_mat(x_coord(i), y_coord(i), :) = est_rate(i,:);
end
%%
cmap = flipud(gray(100)); 
figure; 
colormap(cmap); 
for i = 1:length(est_time)
    image(squeeze(actv_mat(:,:,i))','cdatamapping', 'scaled');
    hold on; 
    scatter(x_coord(viz_pop), y_coord(viz_pop), 10, 'ko');
    hold off;
    set(gca, 'visible', 'off');
    axis square 
    pause(0.2);
    
end

%%
rate_cntr = est_rate(cntpop,:); 
rate_midl = est_rate(midpop,:); 
rate_side = est_rate(sidepop,:);

mean_cntr = mean(rate_cntr,1);
mean_midl = mean(rate_midl,1);
mean_side = mean(rate_side,1);

figure; 
subplot(311); hold on; 
plot(est_time, mean_cntr); 
subplot(312); hold on; 
plot(est_time, mean_midl); 
subplot(313); hold on; 
plot(est_time, mean_side); 

%%
figure; hold on; 
plot(xcorr(mean_cntr, mean_midl, 'coeff'))
plot(xcorr(mean_side, mean_midl, 'coeff'))

%%
x2cntr = round(abs(unq_x - net_arch.centered_neuron.x)) + 1; 
dt_est = est_time(2)-est_time(1);
cmap = parula(max(x2cntr)); 
figure; hold on;
for i = 1:length(x2cntr)
    mean_cntr = mean(est_rate(cntpop,:),1); 
    mean_popi = mean(est_rate(pop_line{i},:),1);
    [corr2cntr, lags] = xcorr(mean_popi, mean_cntr, 'coeff');
    plot(lags*dt_est,corr2cntr, 'Color', cmap(x2cntr(i),:));
end
%%
a = zeros(1,1500);
a(50:end) = 1;
b = zeros(1,1500); 
b(50:end) = 1;

figure; plot(xcorr(a,b,'coeff'))
%% 
[r0,rmax,lags] = xcorr3(actv_mat, est_time);
%%
figure; 
image(lags,'cdatamapping', 'scaled');
axis square