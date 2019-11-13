class nonlinearResistor "Ractual = Rbar + alpha0/(1 + alpha1*Modelica.Math.exp(-alpha2*Vinput)))"
  extends Modelica.Electrical.Analog.Interfaces.OnePort;

  parameter Real alpha0 = 0 "alpha0 non-linear coefficient";
  parameter Real alpha1 = 0 "alpha1 non-linear coefficient";
  parameter Real alpha2 = 0 "alpha2 non-linear coefficient";
  parameter Modelica.SIunits.Resistance Rbar "linear part of the resistance";
  
  final Modelica.SIunits.Resistance Ractual;
  Modelica.Blocks.Interfaces.RealInput Vinput "Pressure that controlled the resistance" annotation(
    Placement(visible = true, transformation(origin = {-92, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 90}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));

  
equation

  Ractual = Rbar + alpha0/(1 + alpha1*Modelica.Math.exp(-alpha2*Vinput));
  v = Ractual * i;
  
  annotation(
    Icon(graphics = {Rectangle(origin = {-6, -31}, lineColor = {0, 0, 255}, extent = {{-54, 59}, {66, 3}}), Line(origin = {-75, 0}, points = {{-15, 0}, {15, 0}}, color = {0, 0, 255}), Line(origin = {-65, 10}, points = {{125, -10}, {155, -10}}, color = {0, 0, 255}), Line(origin = {-65, 10}, points = {{25, -50}, {105, 30}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 7), Text(origin = {0, 65}, lineColor = {0, 0, 255}, extent = {{-60, 15}, {60, -25}}, textString = "%name"), Line(origin = {-65, 10}, points = {{65, 30}, {65, 18}}, color = {0, 0, 255}), Text(origin = {0, -55}, extent = {{-60, 15}, {60, -25}}, textString = "R=%Rbar"), Line(origin = {-55, 20}, points = {{-15, 20}, {55, 20}}, color = {0, 0, 255}), Line(origin = {-45, 30}, points = {{-25, 10}, {-25, 50}}, color = {0, 0, 255})}, coordinateSystem(initialScale = 0.1)));
end nonlinearResistor;