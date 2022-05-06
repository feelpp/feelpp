model BouncingBall
  "The 'classic' bouncing ball model with numerical tolerances"
  parameter Real lambda=2 "Coefficient";
  Real y "y";
initial equation
  y = 1;
equation
  der(y) = -lambda*y;
end BouncingBall;