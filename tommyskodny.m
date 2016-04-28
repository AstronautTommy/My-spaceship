tic
close all
load('KexJobbData');
P = closingPrice((7446:9269),:);

I = find(isnan(P));
for n=1:max(size(I));
    P(I(n)) = P(I(n)-1);
end

dP = [NaN(1,size(P,2)); diff(P)];
DoR = 30; %Days of risk, hur många dagar som ska vara med i beräkning av std.
s=dP.^2;
sequities=zeros(41,13);
sfixedIncome=zeros(41,13);
scommodities=zeros(41,13);
sforeignExc=zeros(41,13);
%MAFAST & MASLOW
for long = 100:5:300;
    for short = 15:2:40;
        %Simpel MACROSS
        smshort=filter(1/short*ones(short,1),1,P);
        smlong=filter(1/long*ones(long,1),1,P);
        sposition= sign(smshort-smlong); %BORDE MAN INTE BEHÖVA TA BORT DET FÖRSTA I DENNA (LIKA LÅNGT SOM "LONG"
        risk=filter(1/DoR*ones(DoR,1),1,s);
%         vinst=sum(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)))); %BLIR RISKEN PÅ FEL DAG NU? RISK ÄR BASERAD PÅ dP OCH INTE P
        medel=mean(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        avik=std(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        sSR=medel./avik*sqrt(252);
        sequities(long/5-19,ceil(short/2)-7) = sum(sSR(1:13));
        sfixedIncome(long/5-19,ceil(short/2)-7) = sum(sSR(14:22));
        scommodities(long/5-19,ceil(short/2)-7) = sum(sSR(23:35));
        sforeignExc(long/5-19,ceil(short/2)-7) = sum(sSR(36:40));
        
        %Linear MACROSS
        lmshort=filter(2*(1:short)/(short*(short+1)),1,P);
        lmlong=filter(2*(1:long)/(long*(long+1)),1,P);
        lposition= sign(lmshort-lmlong); %BORDE MAN INTE BEHÖVA TA BORT DET FÖRSTA I DENNA (LIKA LÅNGT SOM "LONG"
        risk=filter(1/DoR*ones(DoR,1),1,s);
%         vinst=sum(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)))); %BLIR RISKEN PÅ FEL DAG NU? RISK ÄR BASERAD PÅ dP OCH INTE P
        medel=mean(dP(long+1:length(dP),:).*(lposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        avik=std(dP(long+1:length(dP),:).*(lposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        lSR=medel./avik*sqrt(252);
        lequities(long/5-19,ceil(short/2)-7) = sum(lSR(1:13));
        lfixedIncome(long/5-19,ceil(short/2)-7) = sum(lSR(14:22));
        lcommodities(long/5-19,ceil(short/2)-7) = sum(lSR(23:35));
        lforeignExc(long/5-19,ceil(short/2)-7) = sum(lSR(36:40));
        
        %Exponential MACROSS
        emshort=filter(2/(short+1),[1 -(1-2/(short+1))],P); %P(1,:)*(1-2/(short+1)) kan behövas
        emlong=filter(2/(long+1),[1 -(1-2/(long+1))],P); %P(1,:)*(1-2/(long+1)) kan behövas
        eposition= sign(emshort-emlong); %BORDE MAN INTE BEHÖVA TA BORT DET FÖRSTA I DENNA (LIKA LÅNGT SOM "LONG"
        risk=filter(1/DoR*ones(DoR,1),1,s);
%         vinst=sum(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)))); %BLIR RISKEN PÅ FEL DAG NU? RISK ÄR BASERAD PÅ dP OCH INTE P
        medel=mean(dP(long+1:length(dP),:).*(eposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        avik=std(dP(long+1:length(dP),:).*(eposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
        eSR=medel./avik*sqrt(252);
        eequities(long/5-19,ceil(short/2)-7) = sum(eSR(1:13));
        efixedIncome(long/5-19,ceil(short/2)-7) = sum(eSR(14:22));
        ecommodities(long/5-19,ceil(short/2)-7) = sum(eSR(23:35));
        eforeignExc(long/5-19,ceil(short/2)-7) = sum(eSR(36:40));
    end
end
surf(15:2:40, 100:5:300, sequities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('equities with simple MACROSS');
figure;
surf(15:2:40, 100:5:300, lequities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('equities with linear MACROSS');
figure;
surf(15:2:40, 100:5:300, eequities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('equities with exponential MACROSS');
figure;
surf(15:2:40, 100:5:300, sfixedIncome);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('fixedIncome with simple MACROSS');
figure;
surf(15:2:40, 100:5:300, lfixedIncome);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('fixedIncome with linear MACROSS');
figure;
surf(15:2:40, 100:5:300, efixedIncome);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('fixedIncome with exponential MACROSS');
figure;
surf(15:2:40, 100:5:300, scommodities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('commodities with simple MACROSS');
figure;
surf(15:2:40, 100:5:300, lcommodities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('commodities with linear MACROSS');
figure;
surf(15:2:40, 100:5:300, ecommodities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('commodities with exponential MACROSS');
figure;
surf(15:2:40, 100:5:300, sforeignExc);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('foreignExc with simple MACROSS');
figure;
surf(15:2:40, 100:5:300, lforeignExc);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('foreignExc with linear MACROSS');
figure;
surf(15:2:40, 100:5:300, eforeignExc);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('foreignExc with exponential MACROSS');
toc