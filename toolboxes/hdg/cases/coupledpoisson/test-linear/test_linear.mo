within ;
model test_linear
  import Modelica.Constants.pi;
  parameter Modelica.SIunits.Resistance M_Rb = 1;
  parameter Modelica.SIunits.Resistance M_Rout = 1;
  parameter Modelica.SIunits.Capacitance M_Cb = 1;
  parameter Real M_H = 2;  // imposed by the geometry
  parameter Real M_L = 1;  // imposed by the geometry
  parameter Real M_alpha = 10;
  parameter Real M_beta = 0.5;
  parameter Real M_k = 1;
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation (
    Placement(visible = true, transformation(origin={0,-64},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Capacitor Cbuffer(C = M_Cb, v(fixed=true,
        start=0))                                                                            annotation (
    Placement(visible = true, transformation(origin={-40,-18},    extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Resistor Rbuffer(R = M_Rb)  annotation (
    Placement(visible = true, transformation(origin={-80,0},     extent = {{-10, -10}, {10, 10}}, rotation = 0)));

  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_1
    annotation (Placement(transformation(extent={{-40,30},{-20,50}})));
  Modelica.Electrical.Analog.Basic.Resistor R_out(R=M_Rout)
    annotation (Placement(transformation(extent={{-22,-10},{-2,10}})));
  Modelica.Electrical.Analog.Sources.SignalVoltage Pi_out annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-16})));
  Modelica.Blocks.Sources.RealExpression expressionPi_out(y=M_alpha + M_beta*(
        M_H*time + (M_Rb + M_Rout)*M_L*M_L*M_k*time + M_Cb*M_Rout*(M_H + M_Rb*
        M_L*M_L*M_k)))
    annotation (Placement(transformation(extent={{26,28},{48,54}})));
equation
  connect(Cbuffer.n, ground1.p) annotation (
    Line(points={{-40,-28},{-40,-40},{0,-40},{0,-54}},            color = {0, 0, 255}));
  connect(Rbuffer.n, Cbuffer.p) annotation (
    Line(points={{-70,0},{-40,0},{-40,-8}},          color = {0, 0, 255}));
  connect(Pi_1.p, Cbuffer.p)
    annotation (Line(points={{-40,40},{-40,-8}}, color={0,0,255}));
  connect(R_out.p, Cbuffer.p)
    annotation (Line(points={{-22,0},{-40,0},{-40,-8}}, color={0,0,255}));
  connect(R_out.n, Pi_out.p)
    annotation (Line(points={{-2,0},{40,0},{40,-6}}, color={0,0,255}));
  connect(Pi_out.n, ground1.p) annotation (Line(points={{40,-26},{40,-40},{0,-40},
          {0,-54}}, color={0,0,255}));
  connect(expressionPi_out.y, Pi_out.v) annotation (Line(points={{49.1,41},{60,
          41},{60,-16},{52,-16}}, color={0,0,127}));
  annotation (
    uses(Modelica(version="3.2.3")), experiment(
      StopTime=50,
      __Dymola_NumberOfIntervals=1000,
      __Dymola_Algorithm="Dassl"));
end test_linear;
