// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "precAFP3D",
  "ShortName":"precAFP3D",
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
  "Materials":
  {
    "firstMat": // Physical marker
    {
      "name":"firstMat",
      "file":"false",
      "B":"mu_r*mu_0*H:x:y:z:H:mu_0:mu_r"
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
          "expr":"{pi* cos(pi* x)* sin(pi* y)* sin(pi*z),-2*pi* sin(pi* x)* cos(pi* y)* sin(pi*z),pi* sin(pi* x)* sin(pi* y)* cos(pi*z)}:x:y:z"
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

