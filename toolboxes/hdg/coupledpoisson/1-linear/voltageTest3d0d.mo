class voltageTest3d0d 
  extends Modelica.Electrical.Analog.Interfaces.OnePort;

  parameter Modelica.SIunits.Voltage beta;
  parameter Modelica.SIunits.Voltage alpha;
  parameter Real gamma;

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
  final Modelica.SIunits.Voltage Pi2;
  final Modelica.SIunits.Voltage Pi_out;
 
  
equation

  
  Pi1 = beta + alpha*Modelica.Math.sin(gamma*time)+(1-Rb)*Cb*alpha*gamma*Modelica.Math.cos(gamma*time);   
  
  R1 = R1bar + alpha0/(1+alpha1*Modelica.Math.exp(-alpha2*Pi1));
  
  Pi2 = Pi1 + R1*Cb*der(Pi1);

  Pi_out = Pi1 + R1*Cb*der(Pi1) + Rout*Cb *der(Pi1) + Rout*C1*der(Pi2);
  
  v = beta + alpha*Modelica.Math.sin(gamma*time)* ( 1-R1*(1-Rb)*Cb^2*gamma^2 - Rout*(1-Rb)*Cb*gamma^2*(C1+Cb) ) + alpha*gamma*Modelica.Math.cos(gamma*time)*( (1-Rb)*Cb + Rout*C1*(1-R1*(1-Rb)*Cb^2*gamma^2) );
  
annotation(
    Icon(graphics = {Line(origin = {-65, 0}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Ellipse(origin = {-1, 1}, lineColor = {0, 0, 255}, extent = {{-39, 39}, {41, -41}}, endAngle = 360), Line(origin = {65.5029, -0.0767754}, points = {{-25, 0}, {25, 0}, {25, 0}}, color = {0, 0, 255}), Text(origin = {1, -66}, lineColor = {0, 0, 255}, extent = {{-53, 18}, {53, -18}}, textString = "%name"), Text(origin = {-80, 31}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "+", fontSize = 30), Text(origin = {80, 35}, lineColor = {0, 0, 255}, lineThickness = 0.75, extent = {{-44, 17}, {44, -17}}, textString = "-", fontSize = 50), Line(origin = {0, -2}, points = {{20, -22}, {-20, 22}}, color = {0, 0, 255}, arrow = {Arrow.None, Arrow.Filled}, arrowSize = 5)}, coordinateSystem(initialScale = 0.1)));
    
end voltageTest3d0d;
