function pcolor_mike(M,X,Y,ED,CL,DA)
% This function is a custom plot used in Figure 2 of the manuscript


X = [X X(end)+X(2)-X(1)];
Y = [Y Y(end)+Y(2)-Y(1)];

% M = M - min(M(:));
% M = M./max(M(:));
% M = ceil(M*99)+1;
d = 0.07;

for xx = 1:length(X)-1
    for yy = 1:length(Y)-1
        if M(yy,xx)<=DA(1)||M(yy,xx)>=1-DA(1),colorpatch=CL(1,:);else,colorpatch=CL(3,:);end
        pp = patch([X(xx)+d X(xx+1)-d X(xx+1)-d X(xx)+d],Y([yy yy yy+1 yy+1])+[d d -d -d],colorpatch);
        set(pp,'edgecolor','none')
    end
end