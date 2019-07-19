function grad = ExamplecomputeGradient(x)

% Test the "lbfgs" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.


  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);

  grad(1) = -400*x1*(x2-x1^2) - 2*(1-x1);
  grad(2) = 200*(x2-x1^2) + 20.2*(x2-1) + 19.8*(x4-1);
  grad(3) = -360*x3*(x4-x3^2) -2*(1-x3);
  grad(4) = 180*(x4-x3^2) + 20.2*(x4-1) + 19.8*(x2-1);

  