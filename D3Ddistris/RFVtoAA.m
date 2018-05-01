function [n, omega] = RFVtoAA(RFvec)
    rho = norm(RFvec);
    n = RFvec/rho;
    omega = 2*atand(rho);
end