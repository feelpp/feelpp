model test_om
  Modelica.Electrical.Analog.Sources.ConstantVoltage constantVoltage1(V = 5)  annotation(
    Placement(visible = true, transformation(origin = {-18, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Ground ground1 annotation(
    Placement(visible = true, transformation(origin = {24, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Sensors.PotentialSensor potentialSensor1 annotation(
    Placement(visible = true, transformation(origin = {26, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Electrical.Analog.Basic.Resistor resistor1(R = 10)  annotation(
    Placement(visible = true, transformation(origin = {-65, 25}, extent = {{-17, -17}, {17, 17}}, rotation = 180)));
  Modelica.Electrical.Analog.Basic.Ground ground2 annotation(
    Placement(visible = true, transformation(origin = {-6, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(resistor1.n, ground2.p) annotation(
    Line(points = {{-82, 26}, {-96, 26}, {-96, -20}, {-6, -20}, {-6, -20}}, color = {0, 0, 255}));
  connect(constantVoltage1.p, resistor1.p) annotation(
    Line(points = {{-28, 26}, {-48, 26}, {-48, 24}, {-48, 24}, {-48, 26}}, color = {0, 0, 255}));
  connect(potentialSensor1.p, constantVoltage1.p) annotation(
    Line(points = {{16, 56}, {-28, 56}, {-28, 26}}, color = {0, 0, 255}));
  connect(constantVoltage1.n, ground1.p) annotation(
    Line(points = {{-8, 26}, {24, 26}, {24, 24}, {24, 24}}, color = {0, 0, 255}));
  annotation(
    uses(Modelica(version = "3.2.1")));
end test_om;
