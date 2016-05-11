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
risk=filter(1/DoR*ones(DoR,1),1,s);
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
        
        srV=dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)));
%         medel=mean(rV);
%         avik=var(rV);
%         sSR=medel./avik*sqrt(252);
        sett=sum(srV(:,1:13),2);
        stwo=sum(srV(:,14:23),2);
        stre=sum(srV(:,24:35),2);
        sfyr=sum(srV(:,36:40),2);
        sequities(long/5-19,ceil(short/2)-7) = mean(sett)/std(sett)*sqrt(252);
        sfixedIncome(long/5-19,ceil(short/2)-7) = mean(stwo)/std(stwo)*sqrt(252);
        scommodities(long/5-19,ceil(short/2)-7) = mean(stre)/std(stre)*sqrt(252);
        sforeignExc(long/5-19,ceil(short/2)-7) = mean(sfyr)/std(sfyr)*sqrt(252);
        
        %Linear MACROSS
        lmshort=filter(2*(1:short)/(short*(short+1)),1,P);
        lmlong=filter(2*(1:long)/(long*(long+1)),1,P);
        lposition= sign(lmshort-lmlong); %BORDE MAN INTE BEHÖVA TA BORT DET FÖRSTA I DENNA (LIKA LÅNGT SOM "LONG"
%         vinst=sum(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)))); %BLIR RISKEN PÅ FEL DAG NU? RISK ÄR BASERAD PÅ dP OCH INTE P
        lrV=dP(long+1:length(dP),:).*(lposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)));
%         medel=mean(dP(long+1:length(dP),:).*(lposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
%         avik=var(dP(long+1:length(dP),:).*(lposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
%         lSR=medel./avik*sqrt(252);
        lett=sum(lrV(:,1:13),2);
        ltwo=sum(lrV(:,14:23),2);
        ltre=sum(lrV(:,24:35),2);
        lfyr=sum(lrV(:,36:40),2);
        lequities(long/5-19,ceil(short/2)-7) = mean(lett)/std(lett)*sqrt(252);
        lfixedIncome(long/5-19,ceil(short/2)-7) = mean(ltwo)/std(ltwo)*sqrt(252);
        lcommodities(long/5-19,ceil(short/2)-7) = mean(ltre)/std(ltre)*sqrt(252);
        lforeignExc(long/5-19,ceil(short/2)-7) = mean(lfyr)/std(lfyr)*sqrt(252);
        
        %Exponential MACROSS
        emshort=filter(2/(short+1),[1 -(1-2/(short+1))],P); %P(1,:)*(1-2/(short+1)) kan behövas
        emlong=filter(2/(long+1),[1 -(1-2/(long+1))],P); %P(1,:)*(1-2/(long+1)) kan behövas
        eposition= sign(emshort-emlong); %BORDE MAN INTE BEHÖVA TA BORT DET FÖRSTA I DENNA (LIKA LÅNGT SOM "LONG"
%         vinst=sum(dP(long+1:length(dP),:).*(sposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)))); %BLIR RISKEN PÅ FEL DAG NU? RISK ÄR BASERAD PÅ dP OCH INTE P
        erV=dP(long+1:length(dP),:).*(eposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:)));
%         medel=mean(dP(long+1:length(dP),:).*(eposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
%         avik=var(dP(long+1:length(dP),:).*(eposition(long:(length(dP)-1),:)./sqrt(risk(long:length(dP)-1,:))));
%         eSR=medel./avik*sqrt(252);
        eett=sum(erV(:,1:13),2);
        etwo=sum(erV(:,14:23),2);
        etre=sum(erV(:,24:35),2);
        efyr=sum(erV(:,36:40),2);
        eequities(long/5-19,ceil(short/2)-7) = mean(eett)/std(eett)*sqrt(252);
        efixedIncome(long/5-19,ceil(short/2)-7) = mean(etwo)/std(etwo)*sqrt(252);
        ecommodities(long/5-19,ceil(short/2)-7) = mean(etre)/std(etre)*sqrt(252);
        eforeignExc(long/5-19,ceil(short/2)-7) = mean(efyr)/std(efyr)*sqrt(252);
        
        %Bara positiva gissningar, för jämförelse
        
        rV=dP(long+1:length(dP),:)./sqrt(risk(long:length(dP)-1,:));
        ett=sum(rV(:,1:13),2);
        two=sum(rV(:,14:23),2);
        tre=sum(rV(:,24:35),2);
        fyr=sum(rV(:,36:40),2);
        equities(long/5-19,ceil(short/2)-7) = mean(ett)/std(ett)*sqrt(252);
        fixedIncome(long/5-19,ceil(short/2)-7) = mean(two)/std(two)*sqrt(252);
        commodities(long/5-19,ceil(short/2)-7) = mean(tre)/std(tre)*sqrt(252);
        foreignExc(long/5-19,ceil(short/2)-7) = mean(fyr)/std(fyr)*sqrt(252);
    end
end

ind=find(sequities==max(max(sequities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
SEQvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
SEQWIN=sum(SEQvinst,2);

ind=find(sfixedIncome==max(max(sfixedIncome)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
SFIvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
SFIWIN=sum(SFIvinst,2);

ind=find(lcommodities==max(max(lcommodities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LCMvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LCMWIN=sum(LCMvinst,2);

ind=find(lforeignExc==max(max(lforeignExc)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LFEvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LFEWIN=sum(LFEvinst,2);

ind=find(lequities==max(max(lequities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LEQvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LEQWIN=sum(LEQvinst,2);

ind=find(lfixedIncome==max(max(lfixedIncome)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LFIvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LFIWIN=sum(LFIvinst,2);

ind=find(lcommodities==max(max(lcommodities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LCMvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LCMWIN=sum(LCMvinst,2);

ind=find(lforeignExc==max(max(lforeignExc)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
LFEvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
LFEWIN=sum(LFEvinst,2);

ind=find(eequities==max(max(eequities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
EEQvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
EEQWIN=sum(EEQvinst,2);

ind=find(efixedIncome==max(max(efixedIncome)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
EFIvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
EFIWIN=sum(EFIvinst,2);

ind=find(ecommodities==max(max(ecommodities)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
ECMvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
ECMWIN=sum(ECMvinst,2);

ind=find(eforeignExc==max(max(eforeignExc)));
lang=((floor((ind-1)/41)+1)+19)*5;
krt=(ind-floor((ind-1)/41)*41+7)*2;
mshort=filter(1/krt*ones(krt,1),1,P);
mlong=filter(1/lang*ones(lang,1),1,P);
position= sign(mshort-mlong); 
EFEvinst=dP(lang+1:length(dP),:).*(position(lang:(length(dP)-1),:)./sqrt(risk(lang:length(dP)-1,:)));
EFEWIN=sum(EFEvinst,2);

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
figure;
surf(15:2:40, 100:5:300, commodities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('commodities with only positive');
figure;
surf(15:2:40, 100:5:300, foreignExc);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('foreignExc with only positive');
figure;
surf(15:2:40, 100:5:300, fixedIncome);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('fixed income with only positive');
figure;
surf(15:2:40, 100:5:300, equities);
xlabel('Fast moving average [# of days]');
ylabel('Slow moving average [# of days]');
zlabel('Sharpe-ratio');
title('equities with only positive');
toc