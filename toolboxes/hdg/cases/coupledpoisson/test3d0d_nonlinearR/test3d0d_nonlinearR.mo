model test3d0d_nonlinearR
  parameter Modelica.SIunits.Resistance M_Rout = 0.5;
  parameter Modelica.SIunits.Resistance M_Rb = 2;
  parameter Modelica.SIunits.Resistance M_R1 = 0.1;
  parameter Modelica.SIunits.Capacitance M_Cb = 0.5;
  parameter Modelica.SIunits.Capacitance M_C1 = 0.5;
  parameter Real M_alpha0 = 0.0;
  parameter Real M_alpha1 = 0.0;
  parameter Real M_alpha2 = 0.0;
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(visible = true, transformation(origin = {20, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor Rout(R = M_Rout) annotation(
    Placement(visible = true, transformation(origin = {58, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Capacitor Cbuffer(C = M_Cb, v(fixed = true, start = 8.6029)) annotation(
    Placement(visible = true, transformation(origin = {-42, -18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Basic.Resistor Rbuffer(R = M_Rb) annotation(
    Placement(visible = true, transformation(origin = {-60, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_1 annotation(
    Placement(visible = true, transformation(origin = {-32, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Capacitor C(C = M_C1, v(fixed = true, start = 9.2953)) annotation(
    Placement(visible = true, transformation(origin = {20, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_2 annotation(
    Placement(visible = true, transformation(origin = {34, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  nonlinearResistor R(Rbar = M_R1, alpha0 = M_alpha0, alpha1 = M_alpha1, alpha2 = M_alpha2) annotation(
    Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  voltage_nonlinearTest3d0d Pi_out(C1 = M_C1, Cb = M_Cb, R1bar = M_R1, Rb = M_Rb, Rout = M_Rout, alpha = 0.5, alpha0 = M_alpha0, alpha1 = M_alpha1, alpha2 = M_alpha2, beta = 10) annotation(
    Placement(visible = true, transformation(origin = {74, -26}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
equation
  connect(Rbuffer.n, Cbuffer.p) annotation(
    Line(points = {{-50, 10}, {-42, 10}, {-42, -8}, {-42, -8}}, color = {0, 0, 255}));
  connect(Cbuffer.p, R.p) annotation(
    Line(points = {{-42, -8}, {-42, -8}, {-42, 10}, {-20, 10}, {-20, 10}}, color = {0, 0, 255}));
  connect(Pi_1.p, Rbuffer.n) annotation(
    Line(points = {{-42, 50}, {-50, 50}, {-50, 10}, {-50, 10}}, color = {0, 0, 255}));
  connect(Pi_out.n, ground1.p) annotation(
    Line(points = {{74, -36}, {74, -36}, {74, -70}, {20, -70}, {20, -70}}, color = {0, 0, 255}));
  connect(Rout.n, Pi_out.p) annotation(
    Line(points = {{68, 10}, {74, 10}, {74, -16}, {74, -16}}, color = {0, 0, 255}));
  connect(Pi_1.phi, R.Vinput) annotation(
    Line(points = {{-21, 50}, {-18, 50}, {-18, 20}, {-16, 20}}, color = {0, 0, 127}));
  connect(R.n, C.p) annotation(
    Line(points = {{0, 10}, {20, 10}, {20, -10}}, color = {0, 0, 255}));
  connect(C.p, Pi_2.p) annotation(
    Line(points = {{20, -10}, {20, -10}, {20, 50}, {24, 50}, {24, 50}}, color = {102, 102, 102}, pattern = LinePattern.Dash));
  connect(C.p, Rout.p) annotation(
    Line(points = {{20, -10}, {20, -10}, {20, 10}, {48, 10}, {48, 10}}, color = {0, 0, 255}));
  connect(ground1.p, C.n) annotation(
    Line(points = {{20, -70}, {20, -70}, {20, -30}, {20, -30}}, color = {0, 0, 255}));
  connect(Cbuffer.n, ground1.p) annotation(
    Line(points = {{-42, -28}, {-42, -46}, {20, -46}, {20, -70}}, color = {0, 0, 255}));
  annotation(
    uses(Modelica(version = "3.2.2")));
end test3d0d_nonlinearR;
