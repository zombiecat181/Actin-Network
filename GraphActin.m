function GraphActin;
load('Actin Network')
N = size(Warray, 1);
F(N) = struct('cdata',[],'colormap',[]);
for f = 1:N
    W = Warray{f};
    for Q = 1:size(W)
        if Q == 1
           WP = [W(Q,1), W(Q,2), W(Q,4)*cos(W(Q,3))+W(Q,1), W(Q,4)*sin(W(Q,3))+W(Q,2), W(Q,6)];
        else
            WP = [WP; [W(Q,1), W(Q,2), W(Q,4)*cos(W(Q,3))+W(Q,1), W(Q,4)*sin(W(Q,3))+W(Q,2), W(Q,6)]];
        end
    end
    for a = 1:size(WP,1)
        x = [WP(a,1), WP(a,3)];
        y = [WP(a,2), WP(a,4)];
        z = exp(-WP(a, 5)*.05);
        line(x, y, 'Color', [1-z, 1, .4*z]);
        axis([boundaryxminus boundaryxplus boundaryyminus boundaryyplus]);
    end
    drawnow;
    %F(f) = getframe;
    clf
end
save('actin movie');
end