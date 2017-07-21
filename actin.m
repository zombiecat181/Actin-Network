function actin(T)
%define initial variables
%1-xposition, 2-yposition, 3-angle, 4-length, 5-capped?, 6-life, 7-id#, 8-mother 
W = [0 0 14*pi/18 10 0 0 1 1]
%number of actin monomers. range between ~10^5 and 10^6
A = 10000;
%number of branching proteins 
B = 200;
%number of capping proteins
C = 200;
%length of extension of polymer per unit time
v = 0;
%each unit of length is the amount a filament is extended by the binding of
%a single actin monomer or ~4nm
boundaryxminus = -200;
boundaryxplus = 200;
boundaryyminus = -200;
boundaryyplus = 200;
Stats = [A B C v size(W,1) W(1,4) mean(W(1:end,4))];
Warray = {W};
kadd = .5*10^-3;
kbranch = 3*10^-5;
kcap = .0001;
kvan = 0.015;
id = 1;
%loop for each individual time step. Each time step is 
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
        WP = [W(1:N,1), W(1:N,2), diag(W(1:N,4)*(cos(W(1:N,3))).')+W(1:N,1), diag(W(1:N,4)*(sin(W(1:N,3))).')+W(1:N,2)];
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
                W(q,4)=Lsnip;
                W(q,5)= -1*t;
        elseif endpositionx < boundaryyminus
            Lsnip = (boundaryxminus-W(q,1))/cos(W(q,3));
                W = [W;boundaryxplus,Lsnip*sin(W(q,3))+W(q,2),W(q,3),W(q,4)-Lsnip,W(q,5),W(q,6),W(q,7),-1*W(q,7)];
                W(q,4)=Lsnip;
                W(q,5)= -1*t ;   
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
        Nbranch = poissrnd(kbranch*W(q,4)*B);
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
        %check for intersections among filaments
        ProspectiveIntersections = [(W(1:N,1).*tan(W(1:N,3))-W(q,1)*tan(W(q,3))+W(q,2)-W(1:N,2))./(tan(W(1:N,3))-tan(W(q,3)))];
        if WP(q,1) < WP(q,3)
            primero = WP(q,1);
            segundo = WP(q,3);
        else
            primero = WP(q,3);
            segundo = WP(q,1);
        end
        start = WP(1:end,1);
        finish = WP(1:end, 3);
        logico = finish-start;
        logico = logico>0;
        points = [diag(logico*start.')+diag((-1*(logico-1))*finish.'),diag(logico*finish.')+diag((-1*(logico-1))*start.')];
        logicostart = ProspectiveIntersections-points(1:end,1);
        logicostart = logicostart >= 0;
        logicofinish = points(1:end,2)-ProspectiveIntersections;
        logicofinish = logicofinish >=0;
        logico = logicostart.*logicofinish;
        logico(logico==0)=NaN;
        Intersections = ProspectiveIntersections.*logico;
        Intersections = Intersections(Intersections >= primero & Intersections <= segundo);
               % check if filament vanishes
        end
        vanish = rand(size(W,1),1);
        vanish = vanish-(1-exp(-kvan));
        vanish = vanish > 0;
        Aplus = sum(W(1:end,4).*(-1*(vanish-1)));
        Cplusp = W(1:end,5).*(-1*(vanish-1));
        Cplus = sum(Cplusp>0);
        Bplusp = W(1:end,8).*(-1*(vanish-1));
        Bplus = sum(Bplusp>0);
        vanish = vanish.*(linspace(1,size(W,1),size(W,1))).';
        vanish = vanish(vanish>0);
        W = W(vanish,1:end);
    Warray = [Warray; W];
    else
        "the system has collapsed"
        break
    end
    totalL = sum(W(:,4));
    A=A+Aplus-Aminus;
    B=B+Bplus-Bminus;
    C=C+Cplus-Cminus;
    Stats = [Stats; A B C v size(W,1) totalL mean(W(1:end,4))];
end
Warray
Stats
save('Actin Network')
end
            
