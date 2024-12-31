function [d_Prime, response_Bias, F_score, precision, recall, specificity, xi_dp, yi_dp, xi_rb, yi_rb, auc_dp, auc_rb] = aw_SDT(hits, misses, false_alarms, correct_rejections, col_val_line)

% function [d_Prime, response_Bias, F_score, precision, recall, specificity, xi_dp, yi_dp, xi_rb, yi_rb, auc_dp, auc_rb] = aw_SDT(hits, misses, false_alarms, correct_rejections, col_val_line) 
% Calculates the variables and plots of signal detection theory (SDT).
% According to the methods of MX Cohen.
%
% 
% Inputs:
% hits -                              number of hits (correct yes) in binary decision-making paradigm.
% misses -                         number of misses (correct no) in binary decision-making paradigm.
% false_alarms -                number of false alarms (incorrect yes) in binary decision-making paradigm.
% correct_rejections -        number of correct rejections (correct no) in binary decision-making paradigm.
% col_val_line -                 colour of d' and response bias values line. ie. 'k'.
% 
% Outputs:
% d_Prime, response_Bias, F_score, 
% precision, recall, specificity, 
% xi_dp, yi_dp, 
% xi_rb, yi_rb, 
% auc_dp, auc_rb
%

%% calculate measures
% step 1: convert to proportions
p_hit = hits / (hits + misses + eps);
p_fa  =  false_alarms / (false_alarms + correct_rejections + eps);

% step 2: convert proportions to z-score
z_hit = norminv(p_hit);
z_fa  = norminv(p_fa);

% step 3: calculate d'
d_Prime = z_hit - z_fa;

% step 4: calculate response bias
response_Bias = -(z_hit + z_fa) / 2;

% step 5: calculate Fscore
F_score = hits / (hits + (false_alarms + misses) / 2);

% step 6: calculate precision, recall, and specificity
precision = hits / (hits + false_alarms + eps);
recall    = hits / (hits + misses + eps);
specificity = correct_rejections / (correct_rejections + false_alarms + eps);

%% create and plot ROC curves
x_axis  = .01:.01:.99;
dp = norminv(x_axis)' - norminv(x_axis);
rb = -( norminv(x_axis)' + norminv(x_axis) )/2;
colorz = col_val_line;
tol = .01;
col_patch = repmat(.4,[3 1]);
lin_wid = 1.5;

% plot parameters
font_tit = 13; font_ax = font_tit - 2;
mk_sz = 10;

% d': find points and plot isosensitivity curves
idx_dp = find(dp>(d_Prime - tol) & dp<(d_Prime + tol));
[yi_dp,xi_dp] = ind2sub(size(dp),idx_dp);
auc_dp = trapz ([0 x_axis(xi_dp) 1], [0 x_axis(yi_dp) 1]);

% repeat for response bias
idx_rb = find(rb>(response_Bias - tol) & rb<(response_Bias + tol));
[yi_rb,xi_rb] = ind2sub(size(rb),idx_rb);
auc_rb = trapz ([0 x_axis(xi_rb) 1], [0 x_axis(yi_rb) 1]);

% create the 2D spaces
figure('Renderer', 'painters', 'Position', [10 10 600 400]); clf

% plot 1
subplot(121); hold on
% tmp_x = [0 x_axis(xi_dp) 1]; tmp_y = [0 x_axis(yi_dp) 1];
% contourf(x_axis, x_axis, dp, 80,'linecolor','none');
% patch(tmp_x, tmp_y, col_patch,'edgecolor','none');% patch([0 1 1],[0 1 0], col_patch,'EdgeColor','none');
plot([0 x_axis(xi_dp) 1], [0 x_axis(yi_dp) 1],[col_val_line '-' ],'LineWidth',lin_wid);
plot(x_axis(xi_dp), x_axis(yi_dp),[col_val_line '.' ],'markersize',mk_sz,'markerfacecolor',col_val_line);
plot(x_axis,x_axis,'color',col_patch,'LineWidth',lin_wid+.5);
text(.8,.08,['{\itAUC = ' num2str(round(auc_dp,2)) '}'],'horizontalalignment','center');
set(gca, 'xtick',[0 .2 .4 .6 .8 1],'ytick',[0 .2 .4 .6 .8 1]); grid on;
xlabel('{\itp}(FA)','FontSize',font_ax); ylabel('{\itp}(H)','FontSize',font_ax);
title(['{\itd}'' = ' num2str(round(d_Prime,2))],'FontSize',font_tit,'FontWeight','normal');
axis square; box on;
% colormap("gray"); colorbar('southoutside');

subplot(122); hold on
% contourf(x_axis, x_axis, rb, 80,'linecolor','none');
plot([0 x_axis(xi_rb) 1], [1 x_axis(yi_rb) 0],[col_val_line '-' ],'LineWidth',lin_wid);
plot(x_axis(xi_rb), x_axis(yi_rb),[col_val_line '.' ],'markersize',mk_sz,'markerfacecolor',col_val_line);
plot(x_axis,fliplr(x_axis),'color',col_patch,'LineWidth',lin_wid+.5);
text(.2,.08,['{\itAUC = ' num2str(round(auc_rb,2)) '}'],'horizontalalignment','center');
set(gca, 'xtick',[0 .2 .4 .6 .8 1],'ytick',[0 .2 .4 .6 .8 1]); grid on;
xlabel('{\itp}(FA)','FontSize',font_ax); ylabel('{\itp}(H)','FontSize',font_ax);
title(['Response Bias = ' num2str(round(response_Bias,2))],'FontSize',font_tit,'FontWeight','normal');
axis square; box on;
% colorbar('southoutside');

end
