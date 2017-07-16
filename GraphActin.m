function GraphActin;
load('Actin Network')
N = size(Warray, 1);
F(N) = struct('cdata',[],'colormap',[]);
for f = 1:N
    W = Warray{f};
    WP = [W(1:end,1), W(1:end,2), diag(W(1:end,4)*(cos(W(1:end,3))).')+W(1:end,1), diag(W(1:end,4)*(sin(W(1:end,3))).')+W(1:end,2), W(1:end,6)];
        x = [(WP(1:end,1)).';(WP(1:end,3)).'];
        y = [(WP(1:end,2)).';(WP(1:end,4)).'];
        line(x, y, 'Color', [0, 0, 0]);
        axis([-200 200 -200 200]);
    drawnow;
    F(f) = getframe;
    clf
end
save('actin movie');
end