model test3d0d_nonlinearC
  import Modelica.Constants.pi;
  parameter Modelica.SIunits.Resistance M_Rout = 0.5;
  parameter Modelica.SIunits.Resistance M_Rb = 2;
  parameter Modelica.SIunits.Resistance M_R1 = 1 / (4 * pi * pi);
  parameter Modelica.SIunits.Capacitance M_Cb = 1;
  parameter Modelica.SIunits.Capacitance M_C1 = 1;
  parameter Real M_alpha0 = 1/(4*pi*pi);
  parameter Real M_alpha1 = 10.0;
  parameter Real M_alpha2 = 5.0;
  parameter Real M_alpha = 0.5;
  parameter Real M_beta = 10;
  parameter Real M_gamma = 2 * pi;
  parameter Real M_eta = 5;
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(visible = true, transformation(origin = {20, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor Rout(R = M_Rout) annotation(
    Placement(visible = true, transformation(origin = {58, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  voltage_nonlinearC Pi_out(C1bar = M_C1, Cb = M_Cb, R1bar = M_R1, Rb = M_Rb, Rout = M_Rout, alpha = M_alpha, alpha0 = M_alpha0, alpha1 = M_alpha1, alpha2 = M_alpha2, beta = M_beta, eta = M_eta, gamma = M_gamma) annotation(
    Placement(visible = true, transformation(origin = {74, -26}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_1 annotation(
    Placement(visible = true, transformation(origin = {-32, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sensors.PotentialSensor Pi_2 annotation(
    Placement(visible = true, transformation(origin = {38, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor Rbuffer(R = M_Rb) annotation(
    Placement(visible = true, transformation(origin = {-70, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Capacitor Cbuffer(C = M_Cb, v(fixed = true, start = 6.8541)) annotation(
    Placement(visible = true, transformation(origin = {-42, -18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  nonlinearResistor R(Rbar = M_R1, alpha0 = M_alpha0, alpha1 = M_alpha1, alpha2 = M_alpha2) annotation(
    Placement(visible = true, transformation(origin = {-10, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  nonlinearCapacitor C1(Cbar = M_C1, eta = M_eta, v(fixed = true, start = 6.8584)) annotation(
    Placement(visible = true, transformation(origin = {20, -18}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
equation
  connect(Pi_2.p, C1.p) annotation(
    Line(points = {{28, 50}, {20, 50}, {20, -8}}, color = {0, 0, 255}));
  connect(Rout.p, C1.p) annotation(
    Line(points = {{48, 10}, {20, 10}, {20, -8}, {20, -8}}, color = {0, 0, 255}));
  connect(R.n, C1.p) annotation(
    Line(points = {{0, 10}, {20, 10}, {20, -8}, {20, -8}}, color = {0, 0, 255}));
  connect(C1.n, ground1.p) annotation(
    Line(points = {{20, -28}, {20, -28}, {20, -70}, {20, -70}}, color = {0, 0, 255}));
  connect(Pi_1.p, Cbuffer.p) annotation(
    Line(points = {{-42, 50}, {-42, 50}, {-42, -8}, {-42, -8}}, color = {0, 0, 255}));
  connect(Pi_1.phi, R.Vinput) annotation(
    Line(points = {{-21, 50}, {-18, 50}, {-18, 20}, {-16, 20}}, color = {0, 0, 127}));
  connect(Cbuffer.p, R.p) annotation(
    Line(points = {{-42, -8}, {-42, -8}, {-42, 10}, {-20, 10}, {-20, 10}}, color = {0, 0, 255}));
  connect(Rbuffer.n, Cbuffer.p) annotation(
    Line(points = {{-60, 10}, {-42, 10}, {-42, -8}, {-42, -8}}, color = {0, 0, 255}));
  connect(Cbuffer.n, ground1.p) annotation(
    Line(points = {{-42, -28}, {-42, -46}, {20, -46}, {20, -70}}, color = {0, 0, 255}));
  connect(Pi_out.n, ground1.p) annotation(
    Line(points = {{74, -36}, {74, -36}, {74, -70}, {20, -70}, {20, -70}}, color = {0, 0, 255}));
  connect(Rout.n, Pi_out.p) annotation(
    Line(points = {{68, 10}, {74, 10}, {74, -16}, {74, -16}}, color = {0, 0, 255}));
  annotation(
    uses(Modelica(version = "3.2.2")));
end test3d0d_nonlinearC;