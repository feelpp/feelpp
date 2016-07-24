// vim: set ft=javascript:
{
    "Name": "heat",
    "ShortName":"heat",
    "Model":"heat",
    "Parameters":
    {
        "pe":
        {
            "value":"2.0"
        }
    },
    "BoundaryConditions":
    {
        "heat":
        {
            "Neumann":
            {
                "Neuamnn":
                {
                    "expr": "0"
                }
            },
            "Dirichlet":
            {
                "Dirichlet":
                {
                  "expr":"0"
                }
            }
        },
        "fluid":
        {
            "outlet":
            {
                "fluid-outlet":
                {
                    "number":1,                       // number of outlet [default=1]
                    "alemesh_bc":"fixed",             // fixed,free [default=fixed]
                    "type":"windkessel",              // free,windkessel [default=free]
                    "windkessel_coupling":"implicit", // explicit, implicit [default=implicit]
                    "windkessel_Rd":6.2e3,            // resistance distal [default=1.0]
                    "windkessel_Rp":400,              // resistance proximal [default=1.0]
                    "windkessel_Cd":2.72e-4           // capacitance [default=1.0]
                }
            }
        }
    },
    "PostProcess":
    {
        "Fields":["velocity","pressure","displacement"]
    }

}
