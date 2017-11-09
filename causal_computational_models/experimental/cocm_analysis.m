function [ edges, non_edges ] = cocm_analysis( expn, grn )
%does correlation-convergentcrossmap analysis
co = [];
cm = [];
for i=1:size(expn,1)
    for j=1:i-1
        i
        j
        [~, ~ , ~ , SugiY , SugiX , origY , origX] = convergent_cross_map(expn(i,:)',expn(j,:)',1,4);
        SugiX(~isfinite(SugiX)) = 0;
        SugiY(~isfinite(SugiY)) = 0;
        cm1 = corrcoef(SugiX, origX);
        cm2 = corrcoef(SugiY, origY);
        cm = cat(1,cm,[cm1(1,2),i,j,1]);
        cm = cat(1,cm,[cm2(1,2),j,i,1]);
        co1 = abs(xcorr(expn(i,:),expn(j,:),100,'coeff'));
        co1 = co1(:,101:end);
        for l=1:length(co1)
            if isnan(co1(l))
                co1(l) = 0;
            end
        end
        r1 = find(co1==max(co1));
        co2 = abs(xcorr(expn(j,:),expn(i,:),100,'coeff'));
        co2 = co2(:,101:end);
        for l=1:length(co2)
            if isnan(co2(l))
                co2(l) = 0;
            end
        end
        r2 = find(co2==max(co2));
        co = cat(1,co,[max(co1),i,j,r1(1)]);
        co = cat(1,co,[max(co2),j,i,r2(1)]);
    end
end
index = 1;
edges = [];
non_edges = [];
for i=1:size(expn,1)
    for j=1:i-1
        flag1 = 0;
        flag2 = 0;
        for k=1:size(grn,1)
            if all(grn(k,1:2)==[i,j])
                flag1 = 1;
            end
            if all(grn(k,1:2)==[j,i])
                flag2 = 1;
            end
            if flag1 && flag2
                break;
            end
        end
        if flag1
            edges = cat(1,edges,[co(index),cm(index)]);
        else
            non_edges = cat(1,non_edges,[co(index),cm(index)]);
        end
        if flag2
            edges = cat(1,edges,[co(index+1),cm(index+1)]);
        else
            non_edges = cat(1,non_edges,[co(index+1),cm(index+1)]);
        end
        index = index+2;
    end
end
plot(edges(:,1),edges(:,2),'ob');
hold all;
plot(non_edges(:,1),non_edges(:,2),'*r');
end

