function [X_norm, mu, sigma] = featureNormalize(X)

    mu = mean(X);
    X_norm = X - mu;

    sigma = std(X_norm);
    X_norm = X_norm./sigma;

end
