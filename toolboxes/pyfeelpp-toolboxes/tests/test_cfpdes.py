import sys
import feelpp
import feelpp.toolboxes as tb
import feelpp.toolboxes.cfpdes as cfpdes
import pandas as pd

def test_cfpdes():
    e.setConfigFile('cfpdes/square/square2d.cfg')   
    f=cfpdes.simulate(dim=2)

    df = pd.read_csv("cfpdes.heat.measures.csv", sep=",")
    df.head(10)