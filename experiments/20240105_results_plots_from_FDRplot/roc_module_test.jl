function noisy(label; λ=0.0)
    if label
        return 1 - λ * rand()
    else
        return λ * rand()
    end
end

labels = rand(Bool, 10000);

scores(λ) =
    map(labels) do label
        noisy(label, λ=λ)
    end

using ROC

roc_good = roc(scores(0.6), labels, true);
roc_bad = roc(scores(1.0), labels, true);

roc_good
scores(0.6)
labels

area_good = AUC(roc_good)
area_bad = AUC(roc_bad)


using Plots


plot(roc_good, label="good");
plot!(roc_bad, label="bad")

# Try this on roc on the other data....