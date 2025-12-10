function c = analyze_ls(Phi, x)
    % Solve min ||Phi*c - x||_2 stably
    [Q,R] = qr(Phi,0);
    c = R \ (Q' * x(:));
end