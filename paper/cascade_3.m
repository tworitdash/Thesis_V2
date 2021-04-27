function [STT, STR, SRT, SRR] = cascade_3(N, S11, S12, S21, S22, S33, S34, S43, S44, Sl)

    I = eye(length(N), length(N));
    
    U1 = inv(I - S22 * Sl * S33 * Sl);
    U2 = inv(I - S33 * Sl * S22 * Sl);

    STT = (S11 + S12 * Sl * U2 * S33 * Sl * S21);
    STR = (S12 * Sl * U2 * S34);
    SRT = (S43 * Sl * U1 * S21);
    SRR = (S44 + S43 * Sl * U1 * S22 * Sl * S34);

end