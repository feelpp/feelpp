Post-Processing and Visualization
=================================


  # Introduction Introduction

  Once the PDE is solved we usually would like to
  visualize the solution and possibly other data or fields associated
  with the problem. Feel++ provides a very powerful framework for
  post-processing allowing to:
  * visualize time-based data,
  * visualize high order functions,
  * visualize traces,
  * use standard visualization tools such as Paraview, Ensight and Gmsh.

  To achieve this, Feel++ defines a so-called \c Exporter object.

  # VisPrinciples General principles

  The library Feel itself does not have any visualization
  capabilities. However, it provides tools to export both scalar and
  vector fields into two formats: \c EnSight and \c Gmsh. The EnSight format
  can be read by the visualization software EnSight
  http://www.ensight.com or, for instance, by the open source package
  Paraview http://www.paraview.org.

  The choice of format depends on several factors, some of them being
  the robustness/capabilities of the visualization packages associated
  and the type of data to be plotted.

  For first (at most second) degree piecewise polynomials defined in
  straight edge/faces meshes, Paraview (or any other software that
  reads the EnSight format) is a very good choice. A trick to use
  Paraview (or any other visualization software that plots first order
  geometrical and finite elements) to visualize high degree
  polynomials in curved meshes, is to use an interpolation operator
  built on high order nodes associated with the polynomial space,
  see~\cite gpena_cprudhomme_acomen. However, for high degree
  polynomials or meshes with curved elements,
  Gmsh http://geuz.org/gmsh is prefered due to its
  adaptive visualization algorithm for this type of finite/geometrical
  element, see \cite Gmsh.

  To illustrate both approaches in visualizing high degree polynomials
  in curved meshes, we plot, in the Figures below, the nodal
  projection of the function $$f(x,y)=\cos(5x) \sin(5y)$$ in the
  unit circle onto several function spaces. Notice the improved
  ``look`` of the projection using high degree polynomials instead of
  linear projections in a finer mesh. We remark also that the
  difference between the two approaches fades away as we increase the
  degree of the polynomials and the order of the geometrical
  elements. However, the additional cost of building the piecewise
  first order finite element space associated with the finer mesh, the
  calculation of the projection onto this space and the smoothness of
  the graphics makes Gmsh's approach more appealing to visualize this
  type of finite elements. We highlight that these algorithms are
  available for meshes composed only of simplices or quadrilaterals
  (see \ref Notations), in 1D, 2D and 3D.

  <center>
  <table border=0px>
  <tr>
  <td width="15%">\image html circle_p1p1.png "P1" width="15%"</td>
  <td width="15%">\image html circle_p2p2.png</td>
  <td width="15%">\image html circle_p2p2_p1interpolator.png</td>
  </tr>
  <tr>
  <td><center>$$u$$: $$P_1$$, $$P_2/P_2$$ and  $$P_2/P_2$$ with $$P_1$$ interpolation </center></td>
  </tr>
  <tr>
  <td width="15%">\image html circle_p3p3.png</td>
  <td width="15%">\image html circle_p3p3_p1interpolator.png</td>
  <td width="15%">\image html circle_p4p4.png</td>
  <td width="15%">\image html circle_p4p4_p1interpolator.png</td>
  <td width="15%">\image html circle_p5p5.png</td>
  <td width="15%">\image html circle_p5p5_p1interpolator.png</td>
  </tr>
  <tr>
  <td><center>$$u$$: $$P_1$$, $$P_2/P_2$$ and  $$P_2/P_2$$ with $$P_1$$ interpolation </center></td>
  </tr>
  </table>
  </center>


  To further illustrate the capabilities of the Feel's exporter to
  Gmsh, we plot in Figures 1D and 3D functions defined in the $$[-1,1]$$
  interval and the unit sphere, respectively.

  <center>
  <table border=0px>
  <tr>
  <td width="15%">\image html sphere_p1.png</td>
  <td width="15%">\image html sphere_p2.png</td>
  <td width="15%">\image html sphere_p3.png</td>
  <td width="15%">\image html sphere_p4.png</td>
  </tr>
  <tr>
  <td><center>Low to high order visualization of $$u$$</center></td>
  </tr>
  </table>
  </center>

  # Visualization Visualization

  We would like to visualize the function $$u=\sin(\pi x)$$ over
  $$\Omega=\{(x,y) \in \mathbb{R}^2 | x^2 + y^2 < 1\}$$. $$\Omega$$
  is approximated by $$\Omega_h$$.

  To define $$\Omega$$ the code reads
  \snippet myexporter.cpp mesh
and $$u$$ :
  \snippet myexporter.cpp function

  We start with an \c Exporter object that allows to visualize the $$P_1$$ interpolant of $$u$$ over $$\Omega$$.


  <a href="#" class="top">top</a>
  <hr>

  # ExporterReference Reference

  ## ExporterReferenceOptions Options

  \par \c exporter.format
  * Type: multiple choice \c string
  * Values: \c gmsh, \c ensight, \c ensightgold
  * Default value: ensightgold
  * Action: \c exporter.format defines the format to save Feel++ data into.

  \par \c exporter.geometry
  * Type: multiple choice \c string
  * Values: \c change_coords_only, \c change, \c  static
  * Default value: change_coords_only
  * Action: \c exporter.geometry tells to the exporter if the mesh changes over time steps : no
  changes(\c static) coordinates only (\c change_coords_only) or remeshed (\c changes)

  \par \c exporter.fileset
  * Type: \c bool
  * Values: 0, 1
  * Default value: 0
  * Action: \c exporter.fileset=0 save one file per timestep per subdomain,  whereas \c exporter.fileset=1 use one file per subdomain to store all time
  steps. \note This option, \c exporter.fileset=1, reduces tremendously the number of files generated for transient simulations.

  \par \c exporter.prefix
  * Type: \c string
  * Default Value: <empty string>
  * Action: \c exporter.prefix defines the prefix to be user by the exporter. It is especially useful when using multiple exporters and avoid name collision.

  \par \c exporter.directory
  * Type: string
  * Default Value: results
  * Action: \c exporter.directory tells where to export the results to

  \subsubsection ExporterReferenceOptionsEnsightGold Ensight Gold specific options

  \par \c exporter.ensightgold.use-sos
  * Type: \c bool
  * Action: if \c exporter.ensightgold.use-sos=0 multiple case files are handle in first case file else the sos file is used to handle multiple case files

  \par \c exporter.ensightgold.save-face
  * Type: \c bool
  * Action: if \c exporter.ensightgold.save-face=1, the exporter saves mesh and fields on marked faces

