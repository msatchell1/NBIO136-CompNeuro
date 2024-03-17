figure;
for dt = [1e-5, 1e-4, 1e-3]
    tVec = 0: dt : 1;
    randVec = randn(size(tVec))*sqrt(dt);
    randWalk = cumsum(randVec);
    
    hold on;
    plot(tVec, randWalk, DisplayName="dt = "+string(dt))
    ylabel("cumulative sum of noise")
    xlabel("time")
    title("randn*sqrt(dt)")
    legend()
end