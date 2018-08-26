function n = vecnorm(x, p, dim)
    if p ~= 2 
        warning(['This vecnorm script is just for L2 norm!'])
    end
    n = sqrt(sum(x.^2, dim));

end
