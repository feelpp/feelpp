/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim:fenc=utf-8:ft=cpp:et:sw=4:ts=4:sts=4

  This file is part of the Feel library

  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
       Date: 2009-04-03

  Copyright (C) 2009 Université de Grenoble 1

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3.0 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/
/**
   \file parameter.hpp
   \author Christophe Prud'homme <christophe.prudhomme@feelpp.org>
   \date 2009-04-03
 */
#ifndef FEELPP_FEELCORE_PARAMETER_H
#define FEELPP_FEELCORE_PARAMETER_H 1

#include <napp/na.hpp>

namespace Feel
{

namespace na {

using vm = NA::named_argument_t<struct vm_tag>;
using options = NA::named_argument_t<struct options_tag>;
using about = NA::named_argument_t<struct about_tag>;
using prefix = NA::named_argument_t<struct prefix_tag>;
using prefix_overload = NA::named_argument_t<struct prefix_overload_tag>;
using sub = NA::named_argument_t<struct sub_tag>;
using opt = NA::named_argument_t<struct opt_tag>;
using path = NA::named_argument_t<struct path_tag>;
using suffix = NA::named_argument_t<struct suffix_tag>;
using filename = NA::named_argument_t<struct filename_tag>;
using sep = NA::named_argument_t<struct sep_tag>;
using directory = NA::named_argument_t<struct directory_tag>;
using chdir = NA::named_argument_t<struct chdir_tag>;
using config = NA::named_argument_t<struct config_tag>;
using subdir = NA::named_argument_t<struct subdir_tag>;
using format = NA::named_argument_t<struct format_tag>;
using argc = NA::named_argument_t<struct argc_tag>;
using argv = NA::named_argument_t<struct argv_tag>;
using remove = NA::named_argument_t<struct remove_tag>;
using logging = NA::named_argument_t<struct logging_tag>;

using verbose = NA::named_argument_t<struct verbose_tag>;
using threading = NA::named_argument_t<struct threading_tag>;


using matrix = NA::named_argument_t<struct matrix_tag>;
using buildGraphWithTranspose = NA::named_argument_t<struct buildGraphWithTranspose_tag>;
using matrixA = NA::named_argument_t<struct matrixA_tag>;
using matrixB = NA::named_argument_t<struct matrixB_tag>;
using formA = NA::named_argument_t<struct formA_tag>;
using formB = NA::named_argument_t<struct formB_tag>;
using rhs = NA::named_argument_t<struct rhs_tag>;
using solution = NA::named_argument_t<struct solution_tag>;
using prec = NA::named_argument_t<struct prec_tag>;
using transpose = NA::named_argument_t<struct transpose_tag>;
using reuse_prec = NA::named_argument_t<struct reuse_prec_tag>;
using reuse_jac = NA::named_argument_t<struct reuse_jac_tag>;
using maxit = NA::named_argument_t<struct maxit_tag>;
using tolerance = NA::named_argument_t<struct tolerance_tag>;
using rtolerance = NA::named_argument_t<struct rtolerance_tag>;
using atolerance = NA::named_argument_t<struct atolerance_tag>;
using dtolerance = NA::named_argument_t<struct dtolerance_tag>;
using stolerance = NA::named_argument_t<struct stolerance_tag>;
using ksp = NA::named_argument_t<struct ksp_tag>;
using pc = NA::named_argument_t<struct pc_tag>;
using pcfactormatsolverpackage = NA::named_argument_t<struct pcfactormatsolverpackage_tag>;
using constant_null_space = NA::named_argument_t<struct constant_null_space_tag>;
using null_space = NA::named_argument_t<struct null_space_tag>;
using near_null_space = NA::named_argument_t<struct near_null_space_tag>;
using test = NA::named_argument_t<struct test_tag>;
using trial = NA::named_argument_t<struct trial_tag>;
using vector = NA::named_argument_t<struct vector_tag>;
using pattern = NA::named_argument_t<struct pattern_tag>;
using pattern_block = NA::named_argument_t<struct pattern_block_tag>;
using parameter_values = NA::named_argument_t<struct parameter_values_tag>;
using diag_is_nonzero = NA::named_argument_t<struct diag_is_nonzero_tag>;
using block = NA::named_argument_t<struct block_tag>;
using copy_values = NA::named_argument_t<struct copy_values_tag>;
using properties = NA::named_argument_t<struct properties_tag>;
using do_threshold = NA::named_argument_t<struct do_threshold_tag>;
using threshold = NA::named_argument_t<struct threshold_tag>;
using init = NA::named_argument_t<struct init_tag>;
using rowstart = NA::named_argument_t<struct rowstart_tag>;
using colstart = NA::named_argument_t<struct colstart_tag>;
using name = NA::named_argument_t<struct name_tag>;
using label = NA::named_argument_t<struct label_tag>;
using nev = NA::named_argument_t<struct nev_tag>;
using ncv = NA::named_argument_t<struct ncv_tag>;
using mpd = NA::named_argument_t<struct mpd_tag>;
using interval_a = NA::named_argument_t<struct interval_a_tag>;
using interval_b = NA::named_argument_t<struct interval_b_tag>;
using backend = NA::named_argument_t<struct backend_tag>;
using problem = NA::named_argument_t<struct problem_tag>;
using solver = NA::named_argument_t<struct solver_tag>;
using spectrum = NA::named_argument_t<struct spectrum_tag>;
using transform = NA::named_argument_t<struct transform_tag>;
using value_on_diagonal = NA::named_argument_t<struct value_on_diagonal_tag>;
using condense = NA::named_argument_t<struct condense_tag>;
using condenser = NA::named_argument_t<struct condenser_tag>;
using local = NA::named_argument_t<struct local_tag>;
// parameter for exporter
using geo = NA::named_argument_t<struct geo_tag>;
using fileset = NA::named_argument_t<struct fileset_tag>;
// parameter for description of geometries
using h = NA::named_argument_t<struct h_tag>;
using scale = NA::named_argument_t<struct scale_tag>;
using dim = NA::named_argument_t<struct dim_tag>;
using order = NA::named_argument_t<struct order_tag>;
using geo_parameters = NA::named_argument_t<struct geo_parameters_tag>;
using in_memory = NA::named_argument_t<struct in_memory_tag>;
using addmidpoint = NA::named_argument_t<struct addmidpoint_tag>;
using usenames = NA::named_argument_t<struct usenames_tag>;
using xmin = NA::named_argument_t<struct xmin_tag>;
using xmax = NA::named_argument_t<struct xmax_tag>;
using ymin = NA::named_argument_t<struct ymin_tag>;
using ymax = NA::named_argument_t<struct ymax_tag>;
using zmin = NA::named_argument_t<struct zmin_tag>;
using zmax = NA::named_argument_t<struct zmax_tag>;
using nx = NA::named_argument_t<struct nx_tag>;
using ny = NA::named_argument_t<struct ny_tag>;
using nz = NA::named_argument_t<struct nz_tag>;
using refine = NA::named_argument_t<struct refine_tag>;
using update = NA::named_argument_t<struct update_tag>;
using physical_are_elementary_regions = NA::named_argument_t<struct physical_are_elementary_regions_tag>;
using parametricnodes = NA::named_argument_t<struct parametricnodes_tag>;
using force_rebuild = NA::named_argument_t<struct force_rebuild_tag>;
using rebuild = NA::named_argument_t<struct rebuild_tag>;
using shear = NA::named_argument_t<struct shear_tag>;
using recombine = NA::named_argument_t<struct recombine_tag>;
using files_path = NA::named_argument_t<struct files_path_tag>;
using depends = NA::named_argument_t<struct depends_tag>;
using optimize3d_netgen = NA::named_argument_t<struct optimize3d_netgen_tag>;
using pre = NA::named_argument_t<struct pre_tag>;
using post = NA::named_argument_t<struct post_tag>;

// parameter for adapt
using metric = NA::named_argument_t<struct metric_tag>;
using model = NA::named_argument_t<struct model_tag>;
using geotracking = NA::named_argument_t<struct geotracking_tag>;
using snapthickness = NA::named_argument_t<struct snapthickness_tag>;
using statistics = NA::named_argument_t<struct statistics_tag>;
using hmin = NA::named_argument_t<struct hmin_tag>;
using hmax = NA::named_argument_t<struct hmax_tag>;
using collapseOnBoundary = NA::named_argument_t<struct collapseOnBoundary_tag>;
using collapseOnBoundaryTolerance = NA::named_argument_t<struct collapseOnBoundaryTolerance_tag>;
// parameter for xmlParse
using kind = NA::named_argument_t<struct kind_tag>;
using type = NA::named_argument_t<struct type_tag>;
using latex = NA::named_argument_t<struct latex_tag>;
using cmdName = NA::named_argument_t<struct cmdName_tag>;
using values = NA::named_argument_t<struct values_tag>;
using dependencies = NA::named_argument_t<struct dependencies_tag>;
using funcs = NA::named_argument_t<struct funcs_tag>;
using mesh = NA::named_argument_t<struct mesh_tag>;
using geoentity = NA::named_argument_t<struct geoentity_tag>;
using pointset = NA::named_argument_t<struct pointset_tag>;
using desc = NA::named_argument_t<struct desc_tag>;
using desc_lib = NA::named_argument_t<struct desc_lib_tag>;
using shape = NA::named_argument_t<struct shape_tag>;
using convex = NA::named_argument_t<struct convex_tag>;
// project and integrate
using sum = NA::named_argument_t<struct sum_tag>;
using accumulate = NA::named_argument_t<struct accumulate_tag>;
using geomap = NA::named_argument_t<struct geomap_tag>;
using straighten = NA::named_argument_t<struct straighten_tag>;
using expr = NA::named_argument_t<struct expr_tag>;
using grad_expr = NA::named_argument_t<struct grad_expr_tag>;
using div_expr = NA::named_argument_t<struct div_expr_tag>;
using curl_expr = NA::named_argument_t<struct curl_expr_tag>;
using pset = NA::named_argument_t<struct pset_tag>;
using quad = NA::named_argument_t<struct quad_tag>;
using quad1 = NA::named_argument_t<struct quad1_tag>;
using arg = NA::named_argument_t<struct arg_tag>;

using quadptloc = NA::named_argument_t<struct quadptloc_tag>;

using extended_doftable = NA::named_argument_t<struct extended_doftable_tag>;



// orders
using order_u = NA::named_argument_t<struct order_u_tag>;
using order_p = NA::named_argument_t<struct order_p_tag>;

using initial_time = NA::named_argument_t<struct initial_time_tag>;
using final_time = NA::named_argument_t<struct final_time_tag>;
using time_step = NA::named_argument_t<struct time_step_tag>;
using strategy = NA::named_argument_t<struct strategy_tag>;
using steady = NA::named_argument_t<struct steady_tag>;
using restart = NA::named_argument_t<struct restart_tag>;
using reverse = NA::named_argument_t<struct reverse_tag>;
using reverse_load = NA::named_argument_t<struct reverse_load_tag>;
using restart_path = NA::named_argument_t<struct restart_path_tag>;
using restart_at_last_save = NA::named_argument_t<struct restart_at_last_save_tag>;
using rank_proc_in_files_name = NA::named_argument_t<struct rank_proc_in_files_name_tag>;
using freq = NA::named_argument_t<struct freq_tag>;
using n_consecutive_save = NA::named_argument_t<struct n_consecutive_save_tag>;

using markerName = NA::named_argument_t<struct markerName_tag>;
using markerAll = NA::named_argument_t<struct markerAll_tag>;
using marker1 = NA::named_argument_t<struct marker1_tag>;
using marker2 = NA::named_argument_t<struct marker2_tag>;
using marker3 = NA::named_argument_t<struct marker3_tag>;
using marker4 = NA::named_argument_t<struct marker4_tag>;
using marker5 = NA::named_argument_t<struct marker5_tag>;
using marker6 = NA::named_argument_t<struct marker6_tag>;
using marker7 = NA::named_argument_t<struct marker7_tag>;
using marker8 = NA::named_argument_t<struct marker8_tag>;
using marker9 = NA::named_argument_t<struct marker9_tag>;
using marker10 = NA::named_argument_t<struct marker10_tag>;
using marker11 = NA::named_argument_t<struct marker11_tag>;
using marker12 = NA::named_argument_t<struct marker12_tag>;

using domain = NA::named_argument_t<struct domain_tag>;
using image = NA::named_argument_t<struct image_tag>;
using domainSpace = NA::named_argument_t<struct domainSpace_tag>;
using imageSpace = NA::named_argument_t<struct imageSpace_tag>;
using range = NA::named_argument_t<struct range_tag>;
using range_extended = NA::named_argument_t<struct range_extended_tag>;
using element = NA::named_argument_t<struct element_tag>;
using element2 = NA::named_argument_t<struct element2_tag>;
using parameter = NA::named_argument_t<struct parameter_tag>;
using sampling = NA::named_argument_t<struct sampling_tag>;
using context = NA::named_argument_t<struct context_tag>;
using context2 = NA::named_argument_t<struct context2_tag>;
using mpi_communications = NA::named_argument_t<struct mpi_communications_tag>;
using properties_space = NA::named_argument_t<struct properties_space_tag>;

using components = NA::named_argument_t<struct components_tag>;
using periodicity = NA::named_argument_t<struct periodicity_tag>;
using periodic = NA::named_argument_t<struct periodic_tag>;

using collect_garbage = NA::named_argument_t<struct collect_garbage_tag>;

using savehdf5 = NA::named_argument_t<struct savehdf5_tag>;
using partitions = NA::named_argument_t<struct partitions_tag>;
using partition_file = NA::named_argument_t<struct partition_file_tag>;
using respect_partition = NA::named_argument_t<struct respect_partition_tag>;
using rebuild_partitions = NA::named_argument_t<struct rebuild_partitions_tag>;
using rebuild_partitions_filename = NA::named_argument_t<struct rebuild_partitions_filename_tag>;
using worldcomm = NA::named_argument_t<struct worldcomm_tag>;
using worldscomm = NA::named_argument_t<struct worldscomm_tag>;
using parallel = NA::named_argument_t<struct parallel_tag>;
using substructuring = NA::named_argument_t<struct substructuring_tag>;
using structured = NA::named_argument_t<struct structured_tag>;

using jacobian = NA::named_argument_t<struct jacobian_tag>;
using residual = NA::named_argument_t<struct residual_tag>;
using currentElt = NA::named_argument_t<struct currentElt_tag>;
using newElt = NA::named_argument_t<struct newElt_tag>;
using space = NA::named_argument_t<struct space_tag>;
using space2 = NA::named_argument_t<struct space2_tag>;
using initial_theta = NA::named_argument_t<struct initial_theta_tag>;
using min_theta = NA::named_argument_t<struct min_theta_tag>;
using forceRelaxation = NA::named_argument_t<struct forceRelaxation_tag>;

using use_tbb = NA::named_argument_t<struct use_tbb_tag>;
using use_harts = NA::named_argument_t<struct use_harts_tag>;
using grainsize = NA::named_argument_t<struct grainsize_tag>;
using partitioner = NA::named_argument_t<struct partitioner_tag>;

using save = NA::named_argument_t<struct save_tag>;
using ddmethod = NA::named_argument_t<struct ddmethod_tag>;
using penaldir = NA::named_argument_t<struct penaldir_tag>;

using close = NA::named_argument_t<struct close_tag>;

using author = NA::named_argument_t<struct author_tag>;
using task = NA::named_argument_t<struct task_tag>;
using email = NA::named_argument_t<struct email_tag>;
using license = NA::named_argument_t<struct license_tag>;
using copyright = NA::named_argument_t<struct copyright_tag>;
using home = NA::named_argument_t<struct home_tag>;
using bugs = NA::named_argument_t<struct bugs_tag>;
using version = NA::named_argument_t<struct version_tag>;

using points_used = NA::named_argument_t<struct points_used_tag>;
using max_points_used = NA::named_argument_t<struct max_points_used_tag>;
using projection = NA::named_argument_t<struct projection_tag>;

using bc = NA::named_argument_t<struct bc_tag>;
using mu = NA::named_argument_t<struct mu_tag>;
using rho = NA::named_argument_t<struct rho_tag>;
using alpha = NA::named_argument_t<struct alpha_tag>;
using tag = NA::named_argument_t<struct tag_tag>;

// create submesh
using only_on_boundary_faces = NA::named_argument_t<struct only_on_boundary_faces_tag>;
using view = NA::named_argument_t<struct view_tag>;

//using solution = NA::named_argument_t<struct solution_tag>;
using solution_key = NA::named_argument_t<struct solution_key_tag>;
using gradient = NA::named_argument_t<struct gradient_tag>;
using gradient_key = NA::named_argument_t<struct gradient_key_tag>;
using inputs = NA::named_argument_t<struct inputs_tag>;
using script = NA::named_argument_t<struct script_tag>;
using use_script = NA::named_argument_t<struct use_script_tag>;
using compute_pde_coefficients = NA::named_argument_t<struct compute_pde_coefficients_tag>;


using keyword = NA::named_argument_t<struct keyword_tag>;
using repository = NA::named_argument_t<struct repository_tag>;
using physic = NA::named_argument_t<struct physic_tag>;


} // namespace na


inline constexpr auto& _vm = NA::identifier<na::vm>;  // Note: no semicolon
inline constexpr auto& _options = NA::identifier<na::options>;
inline constexpr auto& _about = NA::identifier<na::about>;
inline constexpr auto& _prefix = NA::identifier<na::prefix>;
inline constexpr auto& _prefix_overload = NA::identifier<na::prefix_overload>;
inline constexpr auto& _sub = NA::identifier<na::sub>;
inline constexpr auto& _opt = NA::identifier<na::opt>;
inline constexpr auto& _path = NA::identifier<na::path>;
inline constexpr auto& _suffix = NA::identifier<na::suffix>;
inline constexpr auto& _filename = NA::identifier<na::filename>;
inline constexpr auto& _sep = NA::identifier<na::sep>;
inline constexpr auto& _directory = NA::identifier<na::directory>;
inline constexpr auto& _chdir = NA::identifier<na::chdir>;
inline constexpr auto& _config = NA::identifier<na::config>;
inline constexpr auto& _subdir = NA::identifier<na::subdir>;
inline constexpr auto& _format = NA::identifier<na::format>;
inline constexpr auto& _argc = NA::identifier<na::argc>;
inline constexpr auto& _argv = NA::identifier<na::argv>;
inline constexpr auto& _remove = NA::identifier<na::remove>;
inline constexpr auto& _logging = NA::identifier<na::logging>;

inline constexpr auto& _verbose = NA::identifier<na::verbose>;
inline constexpr auto& _threading = NA::identifier<na::threading>;


inline constexpr auto& _matrix = NA::identifier<na::matrix>;
inline constexpr auto& _buildGraphWithTranspose = NA::identifier<na::buildGraphWithTranspose>;
inline constexpr auto& _matrixA = NA::identifier<na::matrixA>;
inline constexpr auto& _matrixB = NA::identifier<na::matrixB>;
inline constexpr auto& _formA = NA::identifier<na::formA>;
inline constexpr auto& _formB = NA::identifier<na::formB>;
inline constexpr auto& _rhs = NA::identifier<na::rhs>;
inline constexpr auto& _solution = NA::identifier<na::solution>;
inline constexpr auto& _prec = NA::identifier<na::prec>;
inline constexpr auto& _transpose = NA::identifier<na::transpose>;
inline constexpr auto& _reuse_prec = NA::identifier<na::reuse_prec>;
inline constexpr auto& _reuse_jac = NA::identifier<na::reuse_jac>;
inline constexpr auto& _maxit = NA::identifier<na::maxit>;
inline constexpr auto& _tolerance = NA::identifier<na::tolerance>;
inline constexpr auto& _rtolerance = NA::identifier<na::rtolerance>;
inline constexpr auto& _atolerance = NA::identifier<na::atolerance>;
inline constexpr auto& _dtolerance = NA::identifier<na::dtolerance>;
inline constexpr auto& _stolerance = NA::identifier<na::stolerance>;
inline constexpr auto& _ksp = NA::identifier<na::ksp>;
inline constexpr auto& _pc = NA::identifier<na::pc>;
inline constexpr auto& _pcfactormatsolverpackage = NA::identifier<na::pcfactormatsolverpackage>;
inline constexpr auto& _constant_null_space = NA::identifier<na::constant_null_space>;
inline constexpr auto& _null_space = NA::identifier<na::null_space>;
inline constexpr auto& _near_null_space = NA::identifier<na::near_null_space>;
inline constexpr auto& _test = NA::identifier<na::test>;
inline constexpr auto& _trial = NA::identifier<na::trial>;
inline constexpr auto& _vector = NA::identifier<na::vector>;
inline constexpr auto& _pattern = NA::identifier<na::pattern>;
inline constexpr auto& _pattern_block = NA::identifier<na::pattern_block>;
inline constexpr auto& _parameter_values = NA::identifier<na::parameter_values>;
inline constexpr auto& _diag_is_nonzero = NA::identifier<na::diag_is_nonzero>;
inline constexpr auto& _block = NA::identifier<na::block>;
inline constexpr auto& _copy_values = NA::identifier<na::copy_values>;
inline constexpr auto& _properties = NA::identifier<na::properties>;
inline constexpr auto& _do_threshold = NA::identifier<na::do_threshold>;
inline constexpr auto& _threshold = NA::identifier<na::threshold>;
inline constexpr auto& _init = NA::identifier<na::init>;
inline constexpr auto& _rowstart = NA::identifier<na::rowstart>;
inline constexpr auto& _colstart = NA::identifier<na::colstart>;
inline constexpr auto& _name = NA::identifier<na::name>;
inline constexpr auto& _label = NA::identifier<na::label>;
inline constexpr auto& _nev = NA::identifier<na::nev>;
inline constexpr auto& _ncv = NA::identifier<na::ncv>;
inline constexpr auto& _mpd = NA::identifier<na::mpd>;
inline constexpr auto& _interval_a = NA::identifier<na::interval_a>;
inline constexpr auto& _interval_b = NA::identifier<na::interval_b>;
inline constexpr auto& _backend = NA::identifier<na::backend>;
inline constexpr auto& _problem = NA::identifier<na::problem>;
inline constexpr auto& _solver = NA::identifier<na::solver>;
inline constexpr auto& _spectrum = NA::identifier<na::spectrum>;
inline constexpr auto& _transform = NA::identifier<na::transform>;
inline constexpr auto& _value_on_diagonal = NA::identifier<na::value_on_diagonal>;
inline constexpr auto& _condense = NA::identifier<na::condense>;
inline constexpr auto& _condenser = NA::identifier<na::condenser>;
inline constexpr auto& _local = NA::identifier<na::local>;
// parameter for exporter
inline constexpr auto& _geo = NA::identifier<na::geo>;
inline constexpr auto& _fileset = NA::identifier<na::fileset>;
// parameter for description of geometries
inline constexpr auto& _h = NA::identifier<na::h>;
inline constexpr auto& _scale = NA::identifier<na::scale>;
inline constexpr auto& _dim = NA::identifier<na::dim>;
inline constexpr auto& _order = NA::identifier<na::order>;
inline constexpr auto& _geo_parameters = NA::identifier<na::geo_parameters>;
inline constexpr auto& _in_memory = NA::identifier<na::in_memory>;
inline constexpr auto& _addmidpoint = NA::identifier<na::addmidpoint>;
inline constexpr auto& _usenames = NA::identifier<na::usenames>;
inline constexpr auto& _xmin = NA::identifier<na::xmin>;
inline constexpr auto& _xmax = NA::identifier<na::xmax>;
inline constexpr auto& _ymin = NA::identifier<na::ymin>;
inline constexpr auto& _ymax = NA::identifier<na::ymax>;
inline constexpr auto& _zmin = NA::identifier<na::zmin>;
inline constexpr auto& _zmax = NA::identifier<na::zmax>;
inline constexpr auto& _nx = NA::identifier<na::nx>;
inline constexpr auto& _ny = NA::identifier<na::ny>;
inline constexpr auto& _nz = NA::identifier<na::nz>;
inline constexpr auto& _refine = NA::identifier<na::refine>;
inline constexpr auto& _update = NA::identifier<na::update>;
inline constexpr auto& _physical_are_elementary_regions = NA::identifier<na::physical_are_elementary_regions>;
inline constexpr auto& _parametricnodes = NA::identifier<na::parametricnodes>;
inline constexpr auto& _force_rebuild = NA::identifier<na::force_rebuild>;
inline constexpr auto& _rebuild = NA::identifier<na::rebuild>;
inline constexpr auto& _shear = NA::identifier<na::shear>;
inline constexpr auto& _recombine = NA::identifier<na::recombine>;
inline constexpr auto& _files_path = NA::identifier<na::files_path>;
inline constexpr auto& _depends = NA::identifier<na::depends>;
inline constexpr auto& _optimize3d_netgen = NA::identifier<na::optimize3d_netgen>;
inline constexpr auto& _pre = NA::identifier<na::pre>;
inline constexpr auto& _post = NA::identifier<na::post>;

// parameter for adapt
inline constexpr auto& _metric = NA::identifier<na::metric>;
inline constexpr auto& _model = NA::identifier<na::model>;
inline constexpr auto& _geotracking = NA::identifier<na::geotracking>;
inline constexpr auto& _snapthickness = NA::identifier<na::snapthickness>;
inline constexpr auto& _statistics = NA::identifier<na::statistics>;
inline constexpr auto& _hmin = NA::identifier<na::hmin>;
inline constexpr auto& _hmax = NA::identifier<na::hmax>;
inline constexpr auto& _collapseOnBoundary = NA::identifier<na::collapseOnBoundary>;
inline constexpr auto& _collapseOnBoundaryTolerance = NA::identifier<na::collapseOnBoundaryTolerance>;
// parameter for xmlParse
inline constexpr auto& _kind = NA::identifier<na::kind>;
inline constexpr auto& _type = NA::identifier<na::type>;
inline constexpr auto& _latex = NA::identifier<na::latex>;
inline constexpr auto& _cmdName = NA::identifier<na::cmdName>;
inline constexpr auto& _values = NA::identifier<na::values>;
inline constexpr auto& _dependencies = NA::identifier<na::dependencies>;
inline constexpr auto& _funcs = NA::identifier<na::funcs>;
inline constexpr auto& _mesh = NA::identifier<na::mesh>;
inline constexpr auto& _geoentity = NA::identifier<na::geoentity>;
inline constexpr auto& _pointset = NA::identifier<na::pointset>;
inline constexpr auto& _desc = NA::identifier<na::desc>;
inline constexpr auto& _desc_lib = NA::identifier<na::desc_lib>;
inline constexpr auto& _shape = NA::identifier<na::shape>;
inline constexpr auto& _convex = NA::identifier<na::convex>;
// project and integrate
inline constexpr auto& _sum = NA::identifier<na::sum>;
inline constexpr auto& _accumulate = NA::identifier<na::accumulate>;
inline constexpr auto& _geomap = NA::identifier<na::geomap>;
inline constexpr auto& _straighten = NA::identifier<na::straighten>;
inline constexpr auto& _expr = NA::identifier<na::expr>;
inline constexpr auto& _grad_expr = NA::identifier<na::grad_expr>;
inline constexpr auto& _div_expr = NA::identifier<na::div_expr>;
inline constexpr auto& _curl_expr = NA::identifier<na::curl_expr>;
inline constexpr auto& _pset = NA::identifier<na::pset>;
inline constexpr auto& _quad = NA::identifier<na::quad>;
inline constexpr auto& _quad1 = NA::identifier<na::quad1>;
inline constexpr auto& _arg = NA::identifier<na::arg>;

inline constexpr auto& _quadptloc = NA::identifier<na::quadptloc>;

inline constexpr auto& _extended_doftable = NA::identifier<na::extended_doftable>;



// orders
inline constexpr auto& _order_u = NA::identifier<na::order_u>;
inline constexpr auto& _order_p = NA::identifier<na::order_p>;

inline constexpr auto& _initial_time = NA::identifier<na::initial_time>;
inline constexpr auto& _final_time = NA::identifier<na::final_time>;
inline constexpr auto& _time_step = NA::identifier<na::time_step>;
inline constexpr auto& _strategy = NA::identifier<na::strategy>;
inline constexpr auto& _steady = NA::identifier<na::steady>;
inline constexpr auto& _restart = NA::identifier<na::restart>;
inline constexpr auto& _reverse = NA::identifier<na::reverse>;
inline constexpr auto& _reverse_load = NA::identifier<na::reverse_load>;
inline constexpr auto& _restart_path = NA::identifier<na::restart_path>;
inline constexpr auto& _restart_at_last_save = NA::identifier<na::restart_at_last_save>;
inline constexpr auto& _rank_proc_in_files_name = NA::identifier<na::rank_proc_in_files_name>;
inline constexpr auto& _freq = NA::identifier<na::freq>;
inline constexpr auto& _n_consecutive_save = NA::identifier<na::n_consecutive_save>;

inline constexpr auto& _markerName = NA::identifier<na::markerName>;
inline constexpr auto& _markerAll = NA::identifier<na::markerAll>;
inline constexpr auto& _marker1 = NA::identifier<na::marker1>;
inline constexpr auto& _marker2 = NA::identifier<na::marker2>;
inline constexpr auto& _marker3 = NA::identifier<na::marker3>;
inline constexpr auto& _marker4 = NA::identifier<na::marker4>;
inline constexpr auto& _marker5 = NA::identifier<na::marker5>;
inline constexpr auto& _marker6 = NA::identifier<na::marker6>;
inline constexpr auto& _marker7 = NA::identifier<na::marker7>;
inline constexpr auto& _marker8 = NA::identifier<na::marker8>;
inline constexpr auto& _marker9 = NA::identifier<na::marker9>;
inline constexpr auto& _marker10 = NA::identifier<na::marker10>;
inline constexpr auto& _marker11 = NA::identifier<na::marker11>;
inline constexpr auto& _marker12 = NA::identifier<na::marker12>;

inline constexpr auto& _domain = NA::identifier<na::domain>;
inline constexpr auto& _image = NA::identifier<na::image>;
inline constexpr auto& _domainSpace = NA::identifier<na::domainSpace>;
inline constexpr auto& _imageSpace = NA::identifier<na::imageSpace>;
inline constexpr auto& _range = NA::identifier<na::range>;
inline constexpr auto& _range_extended = NA::identifier<na::range_extended>;
inline constexpr auto& _element = NA::identifier<na::element>;
inline constexpr auto& _element2 = NA::identifier<na::element2>;
inline constexpr auto& _parameter = NA::identifier<na::parameter>;
inline constexpr auto& _sampling = NA::identifier<na::sampling>;
inline constexpr auto& _context = NA::identifier<na::context>;
inline constexpr auto& _context2 = NA::identifier<na::context2>;
inline constexpr auto& _mpi_communications = NA::identifier<na::mpi_communications>;
inline constexpr auto& _properties_space = NA::identifier<na::properties_space>;

inline constexpr auto& _components = NA::identifier<na::components>;
inline constexpr auto& _periodicity = NA::identifier<na::periodicity>;
inline constexpr auto& _periodic = NA::identifier<na::periodic>;

inline constexpr auto& _collect_garbage = NA::identifier<na::collect_garbage>;

inline constexpr auto& _savehdf5 = NA::identifier<na::savehdf5>;
inline constexpr auto& _partitions = NA::identifier<na::partitions>;
inline constexpr auto& _partition_file = NA::identifier<na::partition_file>;
inline constexpr auto& _respect_partition = NA::identifier<na::respect_partition>;
inline constexpr auto& _rebuild_partitions = NA::identifier<na::rebuild_partitions>;
inline constexpr auto& _rebuild_partitions_filename = NA::identifier<na::rebuild_partitions_filename>;
inline constexpr auto& _worldcomm = NA::identifier<na::worldcomm>;
inline constexpr auto& _worldscomm = NA::identifier<na::worldscomm>;
inline constexpr auto& _parallel = NA::identifier<na::parallel>;
inline constexpr auto& _substructuring = NA::identifier<na::substructuring>;
inline constexpr auto& _structured = NA::identifier<na::structured>;

inline constexpr auto& _jacobian = NA::identifier<na::jacobian>;
inline constexpr auto& _residual = NA::identifier<na::residual>;
inline constexpr auto& _currentElt = NA::identifier<na::currentElt>;
inline constexpr auto& _newElt = NA::identifier<na::newElt>;
inline constexpr auto& _space = NA::identifier<na::space>;
inline constexpr auto& _space2 = NA::identifier<na::space2>;
inline constexpr auto& _initial_theta = NA::identifier<na::initial_theta>;
inline constexpr auto& _min_theta = NA::identifier<na::min_theta>;
inline constexpr auto& _forceRelaxation = NA::identifier<na::forceRelaxation>;

inline constexpr auto& _use_tbb = NA::identifier<na::use_tbb>;
inline constexpr auto& _use_harts = NA::identifier<na::use_harts>;
inline constexpr auto& _grainsize = NA::identifier<na::grainsize>;
inline constexpr auto& _partitioner = NA::identifier<na::partitioner>;

inline constexpr auto& _save = NA::identifier<na::save>;
inline constexpr auto& _ddmethod = NA::identifier<na::ddmethod>;
inline constexpr auto& _penaldir = NA::identifier<na::penaldir>;

inline constexpr auto& _close = NA::identifier<na::close>;

inline constexpr auto& _author = NA::identifier<na::author>;
inline constexpr auto& _task = NA::identifier<na::task>;
inline constexpr auto& _email = NA::identifier<na::email>;
inline constexpr auto& _license = NA::identifier<na::license>;
inline constexpr auto& _copyright = NA::identifier<na::copyright>;
inline constexpr auto& _home = NA::identifier<na::home>;
inline constexpr auto& _bugs = NA::identifier<na::bugs>;
inline constexpr auto& _version = NA::identifier<na::version>;

inline constexpr auto& _points_used = NA::identifier<na::points_used>;
inline constexpr auto& _max_points_used = NA::identifier<na::max_points_used>;
inline constexpr auto& _projection = NA::identifier<na::projection>;

inline constexpr auto& _bc = NA::identifier<na::bc>;
inline constexpr auto& _mu = NA::identifier<na::mu>;
inline constexpr auto& _rho = NA::identifier<na::rho>;
inline constexpr auto& _alpha = NA::identifier<na::alpha>;
inline constexpr auto& _tag = NA::identifier<na::tag>;

// create submesh
inline constexpr auto& _only_on_boundary_faces = NA::identifier<na::only_on_boundary_faces>;
inline constexpr auto& _view = NA::identifier<na::view>;

//inline constexpr auto& _solution = NA::identifier<na::solution>;
inline constexpr auto& _solution_key = NA::identifier<na::solution_key>;
inline constexpr auto& _gradient = NA::identifier<na::gradient>;
inline constexpr auto& _gradient_key = NA::identifier<na::gradient_key>;
inline constexpr auto& _inputs = NA::identifier<na::inputs>;
inline constexpr auto& _script = NA::identifier<na::script>;
inline constexpr auto& _use_script = NA::identifier<na::use_script>;
inline constexpr auto& _compute_pde_coefficients = NA::identifier<na::compute_pde_coefficients>;

inline constexpr auto& _keyword = NA::identifier<na::keyword>;
inline constexpr auto& _repository = NA::identifier<na::repository>;
inline constexpr auto& _physic = NA::identifier<na::physic>;


} // Feel


#endif /* __feelcore_parameter_H */
