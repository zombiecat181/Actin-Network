function actinv2(T);
%define initial variables
W = [0 0 7*pi/18 10 0 0];
A = 100000;
B = 200;
C = 200;
v = 0;
boundaryxminus = -5000;
boundaryxplus = 5000;
boundaryyminus = -5000;
boundaryyplus = 5000;
Stats = [A B C v size(W,1) 10];
Warray = {W};
kadd = 15*10^-4;
kbranch = 1*10^-6;
kcap = .0001;
kvan = 0.005;
%loop for each individual time step
for t=1:T
    K = [0];
    Aplus = 0;
    Aminus = 0;
    Bplus = 0;
    Bminus = 0;
    Cminus = 0;
    Cplus = 0;
%check if system has collapsed
    N = size(W,1);
    if N >= 1
        for n=1:N
            W(n,6) = W(n,6)+1;
            v = floor(kadd*A);
% add monomers to uncapped filaments
            if W(n,5) == 0 
                W(n,4) = W(n,4)+v;
                Aminus = Aminus + v;
                % check if uncapped filaments should be capped
                P1 = rand; 
                pcap = 1-exp(-kcap*C);
                if P1 <= pcap
                    W(n,5) = 1;
                    Cminus = Cminus+1;
                end
            end 
            % check if filament vanishes
            P2 = rand;
            pvan = 1-exp(-kvan*W(n,6));
             if P2 <= pvan
                    K = [n;K];
                    Bplus = Bplus+1;
                    if W(n,5) ~= 0
                        Cplus = Cplus + 1;
                    end
             end
        %check if filament branches
        P3 = rand;
        pbranch = 1-exp(-kbranch*W(n,4)*B);
        Nbranch = floor(log(P3)/log(pbranch));
        if Nbranch > 0
            Bminus = Bminus+Nbranch;
            for br = 1:Nbranch
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
        end
    for r=1:(size(K)-1)
        Aplus = Aplus+W(K(r),4);
        if size(W,1) == K(r)
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
    totalL = sum(W(:,4));
    A=A+Aplus-Aminus;
    B=B+Bplus-Bminus;
    C=C+Cplus-Cminus;
    Stats = [Stats; A B C v size(W,1) totalL];
end
Warray
Stats
save('Actin Network')
end
            
