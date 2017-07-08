function actin(T)
%define initial variables
%1-xposition, 2-yposition, 3-angle, 4-length, 5-capped?, 6-life, 7-id#, 8-mother 
W = [0 0 14*pi/18 10 0 0 1 1]
A = 1000;
B = 200;
C = 200;
v = 0;
boundaryxminus = -100;
boundaryxplus = 100;
boundaryyminus = -100;
boundaryyplus = 100;
Stats = [A B C v size(W,1) W(1,4)];
Warray = {W};
kadd = 1.2*10^-3;
kbranch = 3*10^-5;
kcap = .0001;
kvan = 0.001;
id = 1
%loop for each individual time step
for t=1:T
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
            q = N+1-n;
            W(q,6) = W(q,6)+1;
            v = kadd*A;
% add monomers to uncapped filaments
            if W(q,5) == 0 
                W(q,4) = W(q,4)+v;
                Aminus = Aminus + v;
                % check if uncapped filaments should be capped
                P1 = rand; 
                pcap = 1-exp(-kcap*C);
                if P1 <= pcap
                    W(q,5) = t;
                    Cminus = Cminus+1;
                end
            end 
             %check if filament exceeds boundary condition
        endpositionx = W(q,1)+W(q,4)*cos(W(q,3));
        endpositiony = W(q,2)+W(q,4)*sin(W(q,3));
        if endpositionx > boundaryxplus
                Lsnip = (boundaryxplus-W(q,1))/cos(W(q,3));
                W = [W;boundaryxminus,Lsnip*sin(W(q,3))+W(q,2),W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6),W(q,7),-1*W(q,7)];
                W(q,4)=Lsnip
                W(q,5)= -1*t
        elseif endpositionx < boundaryyminus
            Lsnip = (boundaryxminus-W(q,1))/cos(W(q,3));
            q
                W = [W;boundaryxplus,Lsnip*sin(W(q,3))+W(q,2),W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6),W(q,7),-1*W(q,7)];
                W(q,4)=Lsnip
                W(q,5)= -1*t    
        end
        endpositiony = W(q,2)+W(q,4)*sin(W(q,3));
        if endpositiony > boundaryyplus
            Lsnip = (boundaryyplus-W(q,2))/sin(W(q,3));
            W = [W; W(q,1)+cos(W(q,3))*Lsnip,boundaryyminus,W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6),W(q,7),-1*W(q,7)];
            W(q,4) = Lsnip;
            W(q,5) = -1*t;
        elseif endpositiony < boundaryyminus
            Lsnip = (boundaryyminus-W(q,2))/sin(W(q,3));
            W = [W; W(q,1)+cos(W(q,3))*Lsnip,boundaryyplus,W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6), W(q,7),-1*W(q,7)];
            W(q,4) = Lsnip;
            W(q,5) = -1*t;
        end
             %check if filament branches
        P3 = rand;
        pbranch = 1-exp(-kbranch*W(q,4)*B);
        Nbranch = floor(log(P3)/log(pbranch));
        if Nbranch > 0 && floor(W(q,4)) > 1
            Bminus = Bminus+Nbranch;
            for br = 1:Nbranch
            L = randi([1,floor(W(q,4))],1);
            x = L*cos(W(q,3))+W(q,1);
            y = L*sin(W(q,3))+W(q,2);
            P4 = rand;
            if P4 <= .5
                theta = W(q,3)+7*pi/18;
            else
                theta = W(q,3)-7*pi/18;
            end
            id = id+1;
            W = [W; x y theta 0 0 0 id W(q,7)];
            end
        end
        % check if filament vanishes
            P2 = rand;
            pvan = 1-exp(-kvan*W(q,6));
             if P2 <= pvan
                    if W(q,8) > 0
                        Bplus = Bplus+1;
                    end
                    Aplus = Aplus+W(q,4);
                    if W(q,5) > 0
                        Cplus = Cplus + 1;
                    end
                    if size(W,1) == q
                    W = W(1:(q-1),1:8);
                    else 
                    W = [W(1:(q-1),1:8);W((q+1):end,1:8)];
                    end
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
            
