// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "Linear 2D",
  "ShortName":"lin_2D",
  "Model":"linear",
  "Description":"Magneto Static",
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
    "COIL": // Physical marker
    {
      "name":"firstMat",
      //"file":"false",
      "B":"mu_0*mu_r*H:x:y:H:mu_r:B:mu_0",
      "mu_r":"1"
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
          "expr":"{0,0,y}:x:y:z"
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
  }
}

