function actinv1(T);
W = [0 0 7*pi/18 10 0 0];
Warray = {W};
pcap = .1;
kbranch = .02;
kvan = .009;
v = 1;
for t=1:T
    Ltot=0;
    N = size(W,1);
    K = 0;
    if N >= 1
        for n=1:N
            W(n,6) = W(n,6)+1;
            if W(n,5) == 0 
                W(n,4) = W(n,4)+v;
                P1 = rand; 
                if P1 <= pcap
                    W(n,5) = 1;
                end
            end 
              P2 = rand;
            pvan = 1-exp(-kvan*W(n,6));
                if P2 <= pvan
                    K = [n;K];
                end
        P3 = rand;
        pbranch = 1-exp(-kbranch*W(n,4));
        if P3 <= pbranch 
            L = randi([1,W(n,4)],1);
            x = L*cos(W(n,3))+W(n,1);
            y = L*sin(W(n,3))+W(n,2);
            P4 = rand;
            if P4 <= .5
                theta = W(n,3)+7*pi/18;
            else
                theta = W(n,3)-7*pi/18;
            end
            W = [W; x y theta 0 0 0];
        end
    end
    for r=1:(size(K)-1)
        if size(W,1) == K(1)
        W = W(1:K(r)-1,1:6);
        else 
        W = [W(1:K(r)-1,1:6);W(K(r)+1:end,1:6)];
        end
    end
    Warray = [Warray; W];
    else
        "the system has collapsed"
        break
    end
end
Warray
save('Actin Network')
end
            
