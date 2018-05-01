O = CrysSym();
objMisO = [80, 3, 1, 0];
g1 = eye(3);
g2 = AAToG(objMisO(1), objMisO(2:4));
[objRFvec, ~] = dgInFZ(g1, g2, O);

[n, omega] = RFVtoAA(objRFvec)
