

function getLegend(n,leg,colors,fontsize,lcn)

% Now Create bar chart with nan so it won't show:
b = bar(nan(n,n));
set(b,{'FaceColor'},colors)
lgd = legend(b,leg,'Location',lcn);
lgd.FontSize=fontsize;
end