within ;
model SimpleBuilding
  "User Guide of the BuildingSystems library Chapter 2: Simple Building"
  extends Modelica.Icons.Example;
  BuildingSystems.Buildings.BuildingTemplates.Building1Zone1DBox building(
    
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.OuterWallSingle2014 constructionWall1,
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.OuterWallSingle2014 constructionWall2,
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.OuterWallSingle2014 constructionWall3,
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.OuterWallSingle2014 constructionWall4,
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.RoofSingle2014 constructionCeiling,
    redeclare BuildingSystems.Buildings.Data.Constructions.Thermal.BasePlateSingle2014 constructionBottom,
    redeclare BuildingSystems.Buildings.Data.Constructions.Transparent.HeatProtectionDoubleGlazingUVal14 constructionWindow1,
    redeclare BuildingSystems.Buildings.Data.Constructions.Transparent.HeatProtectionDoubleGlazingUVal14 constructionWindow2,
    redeclare BuildingSystems.Buildings.Data.Constructions.Transparent.DoubleGlazing constructionWindow3,
    redeclare BuildingSystems.Buildings.Data.Constructions.Transparent.HeatProtectionDoubleGlazingUVal14 constructionWindow4,
    InteriorCeilings= false,
    InteriorWalls=false, calcIdealLoads = false, flexibleOrientation = false, heatSources = false,
    height=3.0,
    heightWindow1=1.0,
    heightWindow2=1.0,
    heightWindow3=1.0,
    heightWindow4=1.0,
    length=9.0, moistureSources = false,
    nZones=1, prescribedAirchange = false, show_TSur = false,width=9.0,
    widthWindow1=3.0,
    widthWindow2=1.0,
    widthWindow3=1.0,
    widthWindow4=1.0)
    annotation (Placement(transformation(extent={{12,-10},{32,10}})));
  BuildingSystems.Buildings.Ambience ambience(
    redeclare block WeatherData = BuildingSystems.Climate.WeatherDataMeteonorm.USA_SanFrancisco_Meteonorm_ASCII, calcLwRad = false,
    nSurfaces=building.nSurfacesAmbience)
    annotation (Placement(transformation(extent={{-30,-10},{-10,10}})));
  Modelica.Blocks.Sources.Constant TSetHeating(
    k=293.15)
    annotation (Placement(transformation(extent={{56,20},{44,32}})));
  Modelica.Blocks.Sources.Constant TSetCooling(
    k=297.15)
    annotation (Placement(transformation(extent={{70,10},{58,22}})));
  Modelica.Blocks.Sources.Constant airchange(
    k=0.5)
    annotation (Placement(transformation(extent={{56,-2},{44,10}})));
equation
  connect(ambience.toSurfacePorts, building.toAmbienceSurfacesPorts)
    annotation (Line(points={{-12,4},{13,4}}, color={0,255,0}));
  connect(ambience.toAirPorts, building.toAmbienceAirPorts)
    annotation (Line(points={{-12,-4},{13,-4}}, color={85,170,255}));
  connect(ambience.TAirRef, building.TAirAmb) annotation (Line(points={{-28.2,7},
          {-32,7},{-32,14},{28.2,14},{28.2,9.8}}, color={0,0,127}));
  connect(ambience.xAir, building.xAirAmb) annotation (Line(points={{-28.2,5},{
          -34,5},{-34,16},{30.4,16},{30.4,9.8}}, color={0,0,127}));
  connect(building.T_setHeating[1], TSetHeating.y) annotation (Line(points={{
          31.8,8},{36,8},{36,26},{40,26},{40,26},{43.4,26}}, color={0,0,127}));
  connect(building.T_setCooling[1], TSetCooling.y) annotation (Line(points={{
          31.8,6},{40,6},{40,16},{57.4,16}}, color={0,0,127}));
  connect(building.airchange[1], airchange.y)
    annotation (Line(points={{31.8,4},{43.4,4},{43.4,4}}, color={0,0,127}));

  annotation (experiment(StartTime=0, StopTime=31536000),
    Icon(coordinateSystem(preserveAspectRatio=false,extent={{-40,-20},{80,40}})),
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-40,-20},{80,40}})),
    uses(Modelica(version="1.12"), BuildingSystems(version="2.0.0-beta")),
  Documentation);
end SimpleBuilding;