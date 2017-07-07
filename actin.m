function actin(T)
%define initial variables
W = [0 0 7*pi/18 10 0 0];
A = 1000;
B = 20;
C = 20;
v = 0;
boundaryxminus = -50;
boundaryxplus = 50;
boundaryyminus = -50;
boundaryyplus = 50;
Stats = [A B C v size(W,1)];
Warray = {W};
kadd = 9*10^-4;
kbranch = 3*10^-5;
kcap = .0001;
kvan = 0.0001;
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
           if endpositiony > boundaryyplus
                Lsnipx = (boundaryxplus-W(q,1))/cos(W(q,3));
                Lsnipy = (boundaryyplus-W(q,2))/sin(W(q,3));
                if Lsnipx > Lsnipy
                    W = [W;W(q,1)+Lsnipy*cos(W(q,3)),boundaryyminus,W(q,3),Lsnipx-Lsnipy,t,W(q,6)];
                    ender = W(end,2)+sin(W(end,3))*W(end,4);
                    W = [W;boundaryxminus,ender,W(q,3),W(q,4)-Lsnipx,W(q,5),W(q,6)];
                    W(q,4) = Lsnipy;
                    W(q,5) = t;
                else
                    W = [W;boundaryxminus,W(q,2)+sin(W(q,3))*Lsnipx,W(q,3),Lsnipy-Lsnipx,t,W(q,6)];
                    ender = W(end,1)+cos(W(end,3))*W(end,4);
                    W = [W;ender,boundaryyminus,W(q,3),W(q,4)-Lsnipy,W(q,5),W(q,6)];
                    W(q,4) = Lsnipx;
                    W(q,5) = t;
                end
            elseif endpositiony < boundaryyminus
                Lsnipx = (boundaryxplus-W(q,1))/cos(W(q,3));
                Lsnipy = (boundaryyminus-W(q,2))/sin(W(q,3));
                if Lsnipx > Lsnipy
                    W = [W;W(q,1)+Lsnipy*cos(W(q,3)),boundaryyplus,W(q,3),Lsnipx-Lsnipy,t,W(q,6)];
                    ender = W(end,2)+sin(W(end,3))*W(end,4);
                    W = [W;boundaryxminus,ender,W(q,3),W(q,4)-Lsnipx,W(q,5),W(q,6)];
                    W(q,4) = Lsnipy;
                    W(q,5) = t;
                else
                    W = [W;boundaryxminus,W(q,2)+Lsnipx*sin(W(q,3)),W(q,3),Lsnipy-Lsnipx,t,W(q,6)];
                    ender = W(end,1)+cos(W(end,3))*W(end,4);
                    W = [W;ender,boundaryyminus,W(q,3),W(q,4)-Lsnipy,W(q,5),W(q,6)];
                    W(q,4) = Lsnipx;
                    W(q,5) = t;
                end
            else
                Lsnip = (boundaryxplus-W(q,1))/cos(W(q,3));
                W = [W;boundaryxminus,Lsnip*sin(W(q,3))+W(q,2),W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6)];
                W(q,4)=Lsnip
                W(q,5)=t
            end
        elseif endpositionx < boundaryxminus
            if endpositiony > boundaryyplus
                Lsnipx = (boundaryxminus-W(q,1))/cos(W(q,3));
                Lsnipy = (boundaryyminus-W(q,2))/sin(W(q,3));
                if Lsnipx > Lsnipy
                    W = [W;W(q,1)+Lsnipy*cos(W(q,3)),boundaryyminus,W(q,3),Lsnipx-Lsnipy,t,W(q,6)];
                    ender = W(end,2)+sin(W(end,3))*W(end,4);
                    W = [W;boundaryxplus,ender,W(q,3),W(q,4)-Lsnipx,W(q,5),W(q,6)];
                    W(q,4) = Lsnipy;
                    W(q,5) = t;
                else
                    W = [W;boundaryxminus,W(q,2)+sin(W(q,3))*Lsnipx,W(q,3),Lsnipy-Lsnipx,t,W(q,6)];
                    ender = W(end,1)+cos(W(end,3))*W(end,4);
                    W = [W;ender,boundaryyminus,W(q,3),W(q,4)-Lsnipy,W(q,5),W(q,6)];
                    W(q,4) = Lsnipx;
                    W(q,5) = t;
                end
            elseif endpositiony < boundaryyminus
                Lsnipx = (boundaryxminus-W(q,1))/cos(W(q,3));
                Lsnipy = (boundaryyminus-W(q,2))/sin(W(q,3));
                if Lsnipx > Lsnipy
                    W = [W;W(q,1)+Lsnipy*cos(W(q,3)),boundaryyplus,W(q,3),Lsnipx-Lsnipy,t,W(q,6)];
                    ender = W(end,2)+sin(W(end,3))*W(end,4);
                    W = [W;boundaryxplus,ender,W(q,3),W(q,4)-Lsnipx,W(q,5),W(q,6)];
                    W(q,4) = Lsnipy;
                    W(q,5) = t;
                else
                    W = [W;boundaryxminus,W(q,2)+sin(W(q,3))*Lsnipx,W(q,3),Lsnipy-Lsnipx,t,W(q,6)];
                    ender = W(end,1)+cos(W(end,3))*W(end,4);
                    W = [W;ender,boundaryyminus,W(q,3),W(q,4)-Lsnipy,W(q,5),W(q,6)];
                    W(q,4) = Lsnipx;
                    W(q,5) = t;
                end
            else
            Lsnip = (-boundaryxminus+W(q,1))/cos(W(q,3));
                W = [W;boundaryxplus,Lsnip*sin(W(q,3))+W(q,2),W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6)];
                W(q,4)=Lsnip
                W(q,5)=t    
            end
        end
        endpositiony = W(q,2)+W(q,4)*sin(W(q,3));
        if endpositiony > boundaryyplus
            Lsnip = (boundaryyplus-W(q,2))/sin(W(q,3));
            W = [W; W(q,1)+cos(W(q,4))*Lsnip,boundaryyminus,W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6)];
            W(q,4) = Lsnip;
            W(q,5) = t;
        elseif endpositiony < boundaryyminus
            Lsnip = (-boundaryyminus+W(q,2))/sin(W(q,3));
            W = [W; W(q,1)+cos(W(q,4))*Lsnip,boundaryyplus,W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6)];
            W(q,4) = Lsnip;
            W(q,5) = t;
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
            W = [W; x y theta 0 0 0];
            end
        end
        % check if filament vanishes
            P2 = rand;
            pvan = 1-exp(-kvan*W(q,6));
             if P2 <= pvan
                    Bplus = Bplus+1;
                    Aplus = Aplus+W(q,4);
                    if W(q,5) > 0
                        Cplus = Cplus + 1;
                    end
                    if size(W,1) == q
                    W = W(1:(q-1),1:6);
                    else 
                    W = [W(1:(q-1),1:6);W((q+1):end,1:6)];
                    end
             end
        end
    Warray = [Warray; W];
    else
        "the system has collapsed"
        break
    end
    A=A+Aplus-Aminus;
    B=B+Bplus-Bminus;
    C=C+Cplus-Cminus;
    Stats = [Stats; A B C v size(W,1)];
end
Warray
Stats
save('Actin Network')
end
            
