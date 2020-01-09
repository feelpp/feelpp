class voltage_nonlinearTest3d0d
  extends Modelica.Electrical.Analog.Interfaces.OnePort;
  import Modelica.Constants.pi;
  parameter Modelica.SIunits.Voltage beta;
  parameter Modelica.SIunits.Voltage alpha;
  Real gamma;

  parameter Modelica.SIunits.Resistance Rb;
  parameter Modelica.SIunits.Resistance R1bar;
  parameter Real alpha0;
  parameter Real alpha1;
  parameter Real alpha2;
  parameter Modelica.SIunits.Resistance Rout;
  parameter Modelica.SIunits.Capacitance Cb;
  parameter Modelica.SIunits.Capacitance C1;
  
  final Modelica.SIunits.Resistance R1;
  final Modelica.SIunits.Voltage Pi1;
  final Modelica.SIunits.Current iR1;
  final Modelica.SIunits.Current i_out;
  // final Modelica.SIunits.Current i_out2;
  final Modelica.SIunits.Voltage Pi2;
  
equation
  gamma = 2*pi;
  Pi1 = beta + alpha*Modelica.Math.sin(gamma*time)+(1-Rb)*Cb*alpha*gamma*Modelica.Math.cos(gamma*time);
  iR1 = (1-Rb)*Cb*Cb*alpha*gamma*gamma*Modelica.Math.sin(gamma*time);
  
  R1 = R1bar + alpha0/(1+alpha1*Modelica.Math.exp(-alpha2*Pi1));

  Pi2 = Pi1 - R1*iR1;
  i_out = C1*der(Pi2)-iR1;
  // i_out2 = C1*der(Pi1) - C1*der(R1)*iR1 -C1*R1*der(iR1) - iR1;
  
  v = Pi2 + Rout * i_out;
  
  annotation(
    Icon(graphics = {Line(origin = {-65, 0}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Ellipse(origin = {-1, 1}, lineColor = {0, 0, 255}, extent = {{-39, 39}, {41, -41}}, endAngle = 360), Line(origin = {65.5029, -0.0767754}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Text(origin = {1, -66}, lineColor = {0, 0, 255}, extent = {{-53, 18}, {53, -18}}, textString = "%name"), Text(origin = {-80, 31}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "+", fontSize = 30), Text(origin = {80, 35}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "-", fontSize = 50), Line(origin = {0, -2}, points = {{20, -22}, {-20, 22}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 5)}, coordinateSystem(initialScale = 0.1)));
end voltage_nonlinearTest3d0d;