function f = ExamplecomputeObjective(x)

% Test the "lbfgs" Matlab interface on the Hock & Schittkowski test problem
% #38. See: Willi Hock and Klaus Schittkowski. (1981) Test Examples for
% Nonlinear Programming Codes. Lecture Notes in Economics and Mathematical
% Systems Vol. 187, Springer-Verlag.

  x1 = x(1);
  x2 = x(2);
  x3 = x(3);
  x4 = x(4);
  
  f = 100*(x2-x1^2)^2 + (1-x1)^2 + 90*(x4-x3^2)^2 + (1-x3)^2 + ...
      10.1*(x2-1)^2 + 10.1*(x4-1)^2 + 19.8*(x2-1)*(x4-1);

