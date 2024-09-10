

region1range = 53:59;
region1length = length(region1range);

region2range = 248:254;
region2length = length(region2range);

concFactor = 10000/300;
% 
% inRegion = sum(kymo(:,248:254),2);
% outRegion = sum(kymo(:,53:59),2);
% figure; plot(inRegion); hold on; plot(outRegion);
%%

colors = [0 0 1; 0 0 0.6; 0 1 0; 0 0.6 0; 1 0 0; 0.6 0 0];
idx = [1 4 2 5 3 6];
figure; hold on;
for i = 1:6

    kymo = allKymos1{idx(i)};
    inRegion = sum(kymo(:,region2range),2);
    plot(inRegion/(region2length*concFactor),'Color',colors(i,:));
end

figure; hold on;
for i = 1:6

    kymo = allKymos1{idx(i)};
    inRegion = sum(kymo(:,region1range),2);
    plot(inRegion/(region1length*concFactor),'Color',colors(i,:));
end


%%

colors = [0 0 1; 0 0 0.6; 0 1 0; 0 0.6 0; 1 0 0; 0.6 0 0];
idx = [1 4 2 5 3 6];
figure; hold on;
for i = 1:6

    kymo = allKymos2{idx(i)};
    inRegion = sum(kymo(:,region1range),2);
    plot(inRegion/(region1length*concFactor),'Color',colors(i,:));
end