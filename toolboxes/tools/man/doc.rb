#!/usr/bin/ruby

#require 'Liquid'

version_major=0
version_minor=107
version_micro=0
version_string="#{version_major}.#{version_minor}"
version_min="v#{version_minor}"

@template = Liquid::Template.parse(File.read(__dir__+"/template.adoc"))

toolboxes =
  { "solid" => "solid",
    "fluid" => "fluid",
    "fsi" => "fsi",
    "hdg_poisson" => "hdg",
    "hdg_coupledpoisson" => "hdg",
    "hdg_elasticity" => "hdg",
    "heat" => "heat",
    "heatfluid" => "heatfluid",
    "electric" => "electric",
    "thermoelectric" => "thermoelectric"}

toolbox_default_cases={ "solid" => "github:{path:toolboxes/solid/cantilever}",
                        "fluid" => "github:{path:toolboxes/fluid/TurekHron}",
                        "fsi" => "github:{path:toolboxes/fsi/TurekHron}",
                        "hdg_poisson" => "github:{path:toolboxes/hdg/poisson/red}",
                        "hdg_coupledpoisson" => "github:{path:toolboxes/hdg/coupledpoisson/1-linear/}",
                        "hdg_elasticity" => "github:{path:toolboxes/hdg/elasticity/quarterturn}",
                        "heat" => "github:{path:toolboxes/heat/Building/ThermalBridgesENISO10211}",
                        "electric" => "github:{path:toolboxes/electric/ElectroMagnets/HL-31_H1}",
                        "heatfluid" => "github:{path:toolboxes/heatfluid/NaturalConvection/cavity}",
                        "thermoelectric" => "github:{path:toolboxes/thermoelectric/ElectroMagnets/HL-31_H1}" }


toolbox_default_cli_cases={ "solid" => "",
                            "fluid" => "",
                            "fsi" => "",
                            "hdg_poisson" => "",
                            "hdg_coupledpoisson" => "",
                            "hdg_elasticity" => "",
                            "heat" => "",
                            "electric" => "",
                            "heatfluid" => "",
                            "thermoelectric" => "" }

toolbox_desc={ "solid" => "solid mechanics",
               "fluid" => "fluid mechanics",
               "fsi" => "fluid structure interaction",
               "hdg_poisson" => "hybridized discontinuous Galerkin Poisson",
               "hdg_coupledpoisson" => "hybridized discontinuous Galerkin Coupled Poisson(3D0D)",
               "hdg_elasticity" => "hybridized discontinuous Galerkin elasticity",
               "heat" => "heat transfer",
               "electric" => "electric",
               "heatfluid" => "heat and fluid",
               "thermoelectric" => "thermoelectric" }

toolbox_docs =
  { "solid" => "csm",
    "fluid" => "cfd",
    "fsi" => "fsi",
    "hdg_poisson" => "",
    "hdg_coupledpoisson" => "",
    "hdg_elasticity" => "",
    "heat" => "heat",
    "heatfluid" => "heatfluid",
    "electric" => "electric",
    "thermoelectric" => "thermoelectric"}

toolbox_default_option_values ={ "solid" => ["3","P1"],
                                 "fluid" => ["3","P2P1G1"],
                                 "fsi" => ["3","P2P1"],
                                 "hdg_poisson" => ["3"],
                                 "hdg_coupledpoisson" => ["3"],
                                 "hdg_elasticity" => ["3"],
                                 "heat" => ["3","P1"],
                                 "electric" => ["3","P1"],
                                 "heatfluid" => ["3","P1-P2P1"],
                                 "thermoelectric" => ["3","P1"] }

toolbox_possible_option_values ={ "solid" => ["2,3","P1,P2"],
                                  "fluid" => ["2,3","P2P1G1,P2P1G2"],
                                  "fsi" => ["2,3","P2P1"],
                                  "hdg_poisson" => ["2,3"],
                                  "hdg_coupledpoisson" => ["2,3"],
                                  "hdg_elasticity" => ["2,3"],
                                  "heat" => ["2,3","P1,P2,P3"],
                                  "electric" => ["2,3","P1,P2,P3"],
                                  "heatfluid" => ["2,3","P1-P2P1"],
                                  "thermoelectric" => ["2,3","P1"] }

toolboxes.each do |app, prefix|
  puts ("toolbox #{app}: write #{prefix}/#{app}.adoc")
  File.delete( "#{prefix}/feelpp_toolbox_#{app}.adoc" ) if File.exist?( "#{prefix}/feelpp_toolbox_#{app}.adoc" )
  File.delete( "#{prefix}/#{app}.adoc" ) if File.exist?( "#{prefix}/#{app}.adoc" )
  File.write("#{prefix}/#{app}.adoc", @template.render(
               "version_min" => "#{version_min}",
               "version_string" => "#{version_string}",
               "toolbox_app" => "feelpp_toolbox_#{app}",
               "toolbox" => "#{app}",
               "toolbox_desc" => toolbox_desc["#{app}"],
               "toolbox_docs" => toolbox_docs["#{app}"],
               "toolbox_default_case" => '"'+toolbox_default_cases["#{app}"]+'"',
               "toolbox_default_cli" => toolbox_default_cli_cases["#{app}"],
               "toolbox_default_option_values" => toolbox_default_option_values["#{app}"],
               "toolbox_possible_option_values" => toolbox_possible_option_values["#{app}"],
             ))
end
