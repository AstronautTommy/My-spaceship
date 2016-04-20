tic
close all
load('KexJobbData');
B = closingPrice((7446:9269),:);

I = find(isnan(B));
for n=1:max(size(I));
    B(I(n)) = B(I(n)-1);
end

% vinststd = zeros(41,13);       %standard deviation on daily returns
% expret = vinststd;                %expected return
% trendr = zeros(41,13); %Correct trend 0/1
% R = vinststd;
% sr = vinststd;                    %sharpe ratio
period=1;
%MAFAST & MASLOW
for int = 0:300:1200;
for q = 100:5:300;
    for b = 15:2:40;
    vinsttyp = zeros(1, 4);
    
    for val = 1:40

%     trendnr = zeros(max(size(T+q:9269-21)),1);
    for n=301+int:601+int
        BMAS(n,val) = mean(B(n-q:n,val));
        BMAF(n,val) = mean(B(n-b:n,val));
        if (ge(BMAF(n,val),BMAS(n,val)))
            a = std(B(n-30:n,val));
            vinsttyp(assetClass(val))= vinsttyp(assetClass(val)) + (B(n+1,val)-B(n,val))./a;
        else
            a = std(B(n-30:n,val));
            vinsttyp(assetClass(val))=vinsttyp(assetClass(val)) + (B(n,val)-B(n+1,val))./a;
        end
%         if (ge(BMAF(n,val),BMAS(n,val))==ge(B(n+21,val)-B(n,val),0));
%             trendnr(n+1-T-q) = 1;
%         end
    end
%     SR(q/5-19,ceil(b/2)-7) = mean(dP)/std(dP);
    end
        equities(q/5-19,ceil(b/2)-7) = vinsttyp(1);
        fixedIncome(q/5-19,ceil(b/2)-7) = vinsttyp(2);
        commodities(q/5-19,ceil(b/2)-7) = vinsttyp(3);
        foreignExc(q/5-19,ceil(b/2)-7) = vinsttyp(4);
    end
end
bestEQ(period)=max(max(equities));
bestFI(period)=max(max(fixedIncome));
bestCM(period)=max(max(commodities));
bestFE(period)=max(max(foreignExc));
period=period+1;
end


surf(15:2:40, 100:5:300, equities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Profit');
title('equities');
figure;
surf(15:2:40, 100:5:300, fixedIncome);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Profit');
title('fixedIncome');
figure;
surf(15:2:40, 100:5:300, commodities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Profit');
title('commodities');
figure;
surf(15:2:40, 100:5:300, foreignExc);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Profit');
title('foreignExc');
toc