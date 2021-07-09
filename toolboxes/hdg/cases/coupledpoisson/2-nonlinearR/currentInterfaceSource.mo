class currentInterfaceSource
  extends Modelica.Electrical.Analog.Interfaces.OnePort;
  import Modelica.Constants.pi;
  parameter Modelica.SIunits.Voltage beta;
  parameter Modelica.SIunits.Voltage alpha;
  parameter Real gamma;

  parameter Modelica.SIunits.Capacitance Cb;
 
  
equation
 
  i = Cb*alpha*gamma*Modelica.Math.cos(gamma*time);
  
  //p.v = beta + alpha*Modelica.Math.sin(gamma*time)+Cb*alpha*gamma*Modelica.Math.cos(gamma*time);
  

  
  annotation(
    Icon(graphics = {Line(origin = {-65, 0}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Ellipse(origin = {-1, 1}, lineColor = {0, 0, 255}, extent = {{-39, 39}, {41, -41}}, endAngle = 360), Line(origin = {65.5029, -0.0767754}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Text(origin = {1, -66}, lineColor = {0, 0, 255}, extent = {{-53, 18}, {53, -18}}, textString = "%name"), Text(origin = {-80, 31}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "+", fontSize = 30), Text(origin = {80, 35}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "-", fontSize = 50), Line(origin = {0, -2}, points = {{20, -22}, {-20, 22}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 5)}, coordinateSystem(initialScale = 0.1)));

end currentInterfaceSource;