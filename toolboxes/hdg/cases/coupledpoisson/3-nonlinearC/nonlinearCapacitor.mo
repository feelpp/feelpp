class nonlinearCapacitor "Cactual = Cbar/(1+eta*Q)"
  extends Modelica.Electrical.Analog.Interfaces.OnePort(v(start=0));

  parameter Real eta "eta non-linear coefficient";
  parameter Modelica.SIunits.Capacitance Cbar "linear part of the capacitance";
  
  Modelica.SIunits.ElectricCharge Q;
  Modelica.SIunits.Capacitance Cactual;

  
equation
  Cactual = Cbar / (1+eta*Q);

  // protect solver from index change
  Q = Cactual*v;
  i = der(Q);
  
  annotation(
    Icon(graphics = {Line(origin = {-75, 0}, points = {{-15, 0}, {55, 0}}, color = {0, 0, 255}), Line(origin = {-65, 10}, points = {{85, -10}, {155, -10}}, color = {0, 0, 255}), Line(origin = {-65, 10}, points = {{25, -50}, {105, 30}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 7), Text(origin = {0, 65}, lineColor = {0, 0, 255}, extent = {{-60, 15}, {60, -25}}, textString = "%name"), Text(origin = {0, -55}, extent = {{-60, 15}, {60, -25}}, textString = "C=%Cbar"), Line(origin = {-20, 1.5}, points = {{0, 40.5}, {0, -41.5}, {0, -39.5}}, color = {0, 0, 255}), Line(origin = {19.9502, 1.05224}, points = {{0, 40.5}, {0, -41.5}, {0, -39.5}}, color = {0, 0, 255})}, coordinateSystem(initialScale = 0.1)));
end nonlinearCapacitor;