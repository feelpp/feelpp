#!/usr/bin/pvpython

import displayGraphics

case="/home/fky/feel/laplacian/dim_1/order_1/beta_1/nu_1/h_0.2/laplacian-1_0.case"
image="./image2D.png"
image2="./image2D2.png"
title="laplacian"
legend="dim=2\norder=1\nh=0.0666667"

displayGraphics.paraviewScreenshot2D(case,image,title,legend,600,1,-45,0,0)
displayGraphics.paraviewScreenshot2D(case,image2,title,legend,600,1,-45,0,0)

case="/home/fky/feel/laplacian/nu_1/beta_1/Simplex_3_1_3/P5/h_0.333333/laplacian-1_0.case"
image="./image3D.png"
title="laplacian"
legend="dim=3\norder=5\nh=0.333333"

#displayGraphics.paraviewScreenshot3D(case,image,title,legend,600,2,-45,0,0)