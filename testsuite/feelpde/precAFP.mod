// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "precAFP",
  "ShortName":"precAFP",
  "Model":"linear",
  "Description":"Magneto Static",
  "Parameters":
  {
    "mu_0":
    {
      "type":"constant",
      "name":"Vacuum permeability",
      //"value":12.5663706144e-7
      "value":1
    }
  },
  "BoundaryConditions":
  {
    "u":
    {
      "Dirichlet":
      {
        "Border":
        {
          "expr":"{0,0}:x:y"
          //"expr":"{pi*sin(pi*x)*cos(pi*y),-pi*cos(pi*x)*sin(pi*y)}:x:y"
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
  } // BoundaryConditions
}

