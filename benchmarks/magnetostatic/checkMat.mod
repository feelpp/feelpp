// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "Team Workshop 13",
  "ShortName":"tws13",
  "Model":"FerroMag",
  "Description":"Ferro Magneto Static",
  "Parameters":
  {
    "mu_0":
    {
      "type":"constant",
      "name":"Vacuum permeability",
      "value":1
    }
  },
  "Materials":
  {
    "Omega":
    {
      "name":"Omega",
      "B": "mu_0*H:x:y:H:mu_0:B",
      "file":"false"
    }
  }, // materials
  "BoundaryConditions":
  {
    "u":
    {
      "Dirichlet":
      {
        "Border":
        {
          "expr":"{0,0,0}:x:y:z"
        }
      }
    },
    "phi":
    {
      "Dirichlet":
      {
        "Border":
        {
          "expr":"0"
        }
      }
    }
  }, // BoundaryConditions
  "PostProcess":
  {
    "Force":["cylinder"],
    "PressureDifference":
    {
      "x1":"{0.15,0.2}",
      "x2":"{0.25,0.2}"
    }
  }
}

