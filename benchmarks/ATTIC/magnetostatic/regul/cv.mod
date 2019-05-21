// -*- mode: javascript -*-  vim: set ft=javascript:
{
  "Name": "Torus Quart",
  "ShortName":"tq",
  "Model":"magnetostatique",
  "Description":"Magnetostatique",
  "Parameters":
  {
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
      "B": "mu_r*mu_0*H:x:y:H:mu_0:mu_r:B",
      "file" : "false",
      "j":"{3*pi*pi*pi* cos(pi*x) *sin(pi *y)* sin(pi* z),-6*pi*pi*pi* sin(pi*x) *cos(pi *y)* sin(pi* z),3*pi*pi*pi* sin(pi*x) *sin(pi *y)* cos(pi* z)}:x:y:z",
      "ex":"{pi* cos(pi* x)* sin(pi* y)* sin(pi*z),-2*pi* sin(pi* x)* cos(pi* y)* sin(pi*z),pi* sin(pi* x)* sin(pi* y)* cos(pi*z)}:x:y:z"
    },
    "Omega":
    {
      "name":"Omega",
      "B": "mu_r*mu_0*H:x:y:H:mu_0:mu_r:B",
      "file":"false",
      "j":"{3*pi*pi*pi* cos(pi*x) *sin(pi *y)* sin(pi* z),-6*pi*pi*pi* sin(pi*x) *cos(pi *y)* sin(pi* z),3*pi*pi*pi* sin(pi*x) *sin(pi *y)* cos(pi* z)}:x:y:z",
      "ex":"{pi* cos(pi* x)* sin(pi* y)* sin(pi*z),-2*pi* sin(pi* x)* cos(pi* y)* sin(pi*z),pi* sin(pi* x)* sin(pi* y)* cos(pi*z)}:x:y:z"
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
    }
  }
}

