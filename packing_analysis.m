
bins = (0:2.5:90);
trialCount = 10;

[~,PSD]=createPorosity_statistics(ones([100,100,100]),xgrid,ygrid,zgrid,Tm,maxTemps,30,bins,trialCount);

figure
for i = 1:trialCount
    plot(bins(2:end)*2,PSD(:,i),"LineWidth",2)
    hold on
end
xlabel("Particle Diameter")
ylabel("Probability")

%%%

writematrix(horzcat(bins(2:end)'*2,PSD),"powder_PSD.csv")

%%%

figure
PSD_vol = PSD;
for i = 1:trialCount
    totalSize = PSD(:,i).*bins(2:end).^3;
    totalVol = sum(totalSize);
    plot(bins(2:end)*2,totalSize/totalVol,"LineWidth",2)
    PSD_vol(:,i)= totalSize/totalVol;
    hold on
end
xlabel("Particle Diameter")
ylabel("Volume Fraction")

writematrix(horzcat(bins(2:end)'*2,PSD_vol*100),"powder_PSD_volumeFrac.csv")

