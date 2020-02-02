#!/usr/bin/ruby

#require 'Liquid'

version_major=0
version_minor=107
version_micro=0
version_string="v#{version_major}.#{version_minor}.#{version_micro}"
version_min="v#{version_minor}"

@template = Liquid::Template.parse(File.read(__dir__+"/template.adoc"))

toolboxes = [ "solid", "fluid", "fsi", "hdg", "heat", "heatfluid", "electric", "thermoelectric" ]

toolbox_default_cases={ "solid" => "github:{path:toolboxes/solid/cantilever}",
                        "fluid" => "github:{path:toolboxes/fluid/TurekHron}",
                        "fsi" => "github:{path:toolboxes/fsi/TurekHron}",
                        "hdg" => "github:{path:toolboxes/hdg/poisson/red}",
                        "heat" => "github:{path:toolboxes/heat/Building/ThermalBridgesENISO10211}",
                        "electric" => "github:{path:toolboxes/electric/ElectroMagnets/HL-31_H1}",
                        "heatfluid" => "github:{path:toolboxes/heatfluid/NaturalConvection/cavity}",
                        "thermoelectric" => "github:{path:toolboxes/thermoelectric/ElectroMagnets/HL-31_H1}" }


toolbox_default_cli_cases={ "solid" => "",
                            "fluid" => "",
                            "fsi" => "",
                            "hdg" => "",
                            "heat" => "",
                            "electric" => "",
                            "heatfluid" => "",
                            "thermoelectric" => "" }

toolbox_desc={ "solid" => "solid mechanics",
               "fluid" => "fluid mechanics",
               "fsi" => "fluid structure interaction",
               "hdg" => "hybridized discontinuous Galerkin",
               "heat" => "heat transfer",
               "electric" => "electric",
               "heatfluid" => "heat and fluid",
               "thermoelectric" => "thermoelectric" }

toolbox_default_option_values ={ "solid" => ["3",""],
                                 "fluid" => ["3",""],
                                 "fsi" => ["3",""],
                                 "hdg" => ["3"],
                                 "heat" => ["3","P1"],
                                 "electric" => ["3","P1"],
                                 "heatfluid" => ["3","P1"],
                                 "thermoelectric" => ["3","P1"] }

toolbox_possible_option_values ={ "solid" => ["2,3",""],
                                  "fluid" => ["2,3",""],
                                  "fsi" => ["2,3",""],
                                  "hdg" => ["2,3"],
                                  "heat" => ["2,3","P1"],
                                  "electric" => ["2,3","P1"],
                                  "heatfluid" => ["2,3","P1"],
                                  "thermoelectric" => ["2,3","P1"] }

for l in toolboxes
  #puts l
  puts ("toolbox #{l}: write #{l}/#{l}.adoc")
  File.delete( "#{l}/feelpp_toolbox_#{l}.adoc" ) if File.exist?( "#{l}/feelpp_toolbox_#{l}.adoc" )
  File.delete( "#{l}/#{l}.adoc" ) if File.exist?( "#{l}/#{l}.adoc" )
  File.write("#{l}/#{l}.adoc", @template.render(
               "version_min" => "#{version_min}",
               "toolbox_app" => "feelpp_toolbox_#{l}",
               "toolbox" => "#{l}",
               "toolbox_desc" => toolbox_desc["#{l}"],
               "toolbox_default_case" => '"'+toolbox_default_cases["#{l}"]+'"',
               "toolbox_default_cli" => toolbox_default_cli_cases["#{l}"],
               "toolbox_default_option_values" => toolbox_default_option_values["#{l}"],
               "toolbox_possible_option_values" => toolbox_possible_option_values["#{l}"],
             ))
end
