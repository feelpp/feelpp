model test3d0d
  import Modelica.Constants.pi;
  parameter Modelica.SIunits.Resistance M_Rb = 1;
  parameter Modelica.SIunits.Capacitance M_Cb = 1;
  parameter Real M_H = 2;
  parameter Real M_L = 1;
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation (
    Placement(visible = true, transformation(origin={0,-64},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Capacitor Cbuffer(C = M_Cb, v(fixed=false,
        start=10))                                                                           annotation (
    Placement(visible = true, transformation(origin={-40,-18},    extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Resistor Rbuffer(R = M_Rb)  annotation (
    Placement(visible = true, transformation(origin={-80,0},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage(V=M_H + M_L*M_L)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}})));

  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_1
    annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
equation
  connect(Cbuffer.n, ground1.p) annotation (
    Line(points={{-40,-28},{-40,-40},{0,-40},{0,-54}},            color = {0, 0, 255}));
  connect(Rbuffer.n, Cbuffer.p) annotation (
    Line(points={{-70,0},{-40,0},{-40,-8}},          color = {0, 0, 255}));
  connect(constantVoltage.p, Cbuffer.p)
    annotation (Line(points={{-10,0},{-40,0},{-40,-8}}, color={0,0,255}));
  connect(constantVoltage.n, ground1.p) annotation (Line(points={{10,0},{20,0},{
          20,-40},{0,-40},{0,-54}}, color={0,0,255}));
  connect(Pi_1.p, Cbuffer.p)
    annotation (Line(points={{-40,40},{-40,-8}}, color={0,0,255}));
  annotation (
    uses(Modelica(version="3.2.3")));
end test3d0d;