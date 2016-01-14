// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "Torus Quart",
  "ShortName":"tq",
  "Model":"magnetostatique",
  "Description":"Magnetostatique",
  "Parameters":
  {
    "mu_0":
    {
      "type":"constant",
      "name":"model parameter 1",
      "value":12.566370614e-7
    },
    "a":
    {
      "type":"constant",
      "name":"model parameter 1",
      "value":1
    },
    "b":
    {
      "type":"constant",
      "name":"model parameter 2",
      "value":5e-3
    },
    "c":
    {
      "type":"constant",
      "name":"model parameter 3",
      "value":5e-3
    }
  },
  "Materials":
  {
    "COIL": // Physical marker
    {
      "name":"COIL", // Name of material
      "j":"{-48.e+6*(0.5/(2*Pi))*y/(x^2+y^2),48.e+6*(0.5/(2*Pi))*x/(x^2+y^2),0}:x:y:z",
      "ex":"{1,1,1}:x:y:z"
    },
    "Omega":
    {
      "name":"Omega",
      "j":"{0,0,0}:x:y:z",
      "ex":"{1,1,1}:x:y:z"
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
          "expr":"0:x:y:z"
        }
      }
    }
  }
}

