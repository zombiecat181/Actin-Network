function actinstats; 
load ('Actin Network');
x = 1:size(Stats);
y1 = Stats(1:end, 1);
y2 = Stats(1:end, 2);
y3 = Stats(1:end, 3);
y4 = Stats(1:end, 4);
y5 = Stats(1:end, 5);
y6 = Stats(1:end, 6);
figure 
subplot(5,1,1)
plot(x,y1)
title('Actin Monomers in Solution')

subplot(5,1,2)
plot(x,y2)
title('Branch Protein in Solution')

subplot(5,1,3)
plot(x,y3)
title('Capping Protein in Solution')

subplot(5,1,4)
plot(x,y4)
title('Velocity')

subplot(5,1,5)
plot(x,y6)
title('total length')

