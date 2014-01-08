# $id: $

# setting like in dune-web
set(CURDIR doc/pdelab-howto)
set(BASEDIR ../..)

# install the html pages
set(docdir ${CMAKE_INSTALL_DOCDIR})
if( BUILD_DOCS)
  set(DOCFILES  ${CMAKE_CURRENT_BINARY_DIR}/pdelab-howto-article.pdf) # pdelab-howto-article.pdf 
  set(EXTRA_DIST  ${DOCFILES})
  set(EXTRAINSTALL  ${DOCFILES})
endif( BUILD_DOCS )

if( BUILD_DOCS)
if(CONVERT_FOUND)
 set(DOCS  ${HTDOCS} ${MAINDOCS} ${WMLDOCS})
endif(CONVERT_FOUND )
endif( BUILD_DOCS )

set(TEX_FILES adaptivity.tex algebraicformulation.tex constrainedspaces.tex 
	functionspaces.tex istl-backend.tex pdelab-howto-beamer.tex 
        summary.tex instationaryproblems.tex introduction.tex istl-backend.tex
	operators.tex parallel_computing.tex
        pdelab-howto-main.tex systems.tex introduction.tex 
        pdelab-howto-article.tex  stationaryproblems.tex summary.tex
	systems.tex weightedresidual.tex workflow.tex)

set(EPS_IMAGES
  EPS/adaptivity/adaptive_iso.eps
  EPS/adaptivity/adaptivity1.eps
  EPS/adaptivity/adaptivity2.eps
  EPS/adaptivity/adaptivity3.eps
  EPS/adaptivity/adaptivity4.eps
  EPS/adaptivity/adaptivity5.eps
  EPS/adaptivity/bisection.eps
  EPS/adaptivity/cartesian.eps
  EPS/adaptivity/errind_maxconc.eps
  EPS/adaptivity/fine_iso.eps
  EPS/adaptivity/hadapt1.eps
  EPS/adaptivity/hadapt.eps
  EPS/adaptivity/hanging.eps
  EPS/adaptivity/hexclosurenew.eps
  EPS/adaptivity/high_resolution.eps
  EPS/adaptivity/low_resolution.eps
  EPS/adaptivity/red_green.eps
  EPS/adaptivity/redgreen.eps
  EPS/adaptivity/refinehex.eps
  EPS/adaptivity/refinepri.eps
  EPS/adaptivity/refinepyr.eps
  EPS/adaptivity/refinetet.eps 
  EPS/aggregates_problem_f_2d_65536_cl015_var8_mean0.eps
  EPS/backend.eps
  EPS/blockstructure.eps
  EPS/crank/crank_DirBC.eps
  EPS/crank/crank_groups.eps
  EPS/crank/crank_sol.eps
  EPS/crank/crank_surfgroups1.eps
  EPS/crank/crank_surfgroups2.eps
  EPS/crank/geo1.eps
  EPS/crank/geo2.eps
  EPS/dd.eps
  EPS/dnapl-3d-het-iso.eps
  EPS/example01a_Q1.eps
  EPS/example01a_Q2.eps
  EPS/example01a_RT.eps
  EPS/example02_Q1.eps
  EPS/example03.eps
  EPS/example04.eps
  EPS/example05_0_05.eps
  EPS/example05_0_50.eps
  EPS/example05_0_5.eps
  EPS/example05_1_05.eps
  EPS/example05_1_50.eps
  EPS/example05_1_5.eps
  EPS/function3.eps
  EPS/function3_fm.eps
  EPS/gcad3d_cyl_dune.eps
  EPS/gcad3d_cyl_gmsh.eps
  EPS/gcad3d_deckel.eps
  EPS/hammer_surfgroups.eps
  EPS/istlsolver.eps
  EPS/iwrlogo.eps
  EPS/iwrlogo-klein.eps
  EPS/konstr_n_uberl_str.eps
  EPS/konstr_n_uberl_unstr.eps
  EPS/konstr_uberl_str.eps
  EPS/konstr_uberl_unstr.eps
  EPS/L/L_geoslayout.eps
  EPS/L/L_geostatistics.eps
  EPS/L/L_internalfaces.eps
  EPS/L/L_mesh.eps
  EPS/L/L_meshstatistics.eps
  EPS/L/L_numbers_vols_faces.eps
  EPS/L/L_physgroups.eps
  EPS/L/L_surfgroup1.eps
  EPS/L/L_surfgroup2.eps
  EPS/L/L_volgroup1.eps
  EPS/L/L_volgroup2.eps
  EPS/modules.eps
  EPS/motor.eps
  EPS/numalayout.eps
  EPS/P1P1b.eps
  EPS/P1P1.eps
  EPS/P2P2b.eps
  EPS/P2P2.eps
  EPS/paket.eps
  EPS/partitionsingle.eps
  EPS/pointblockmatrix.eps
  EPS/q1constrainedinterpolate.eps
  EPS/q1.eps
  EPS/q1interpolate.eps
  EPS/q1laplace3d.eps
  EPS/q1laplace.eps
  EPS/quadrilateral.eps
  EPS/rotbilinear.eps
  EPS/scalarmatrix.eps
  EPS/SecMesh2D.eps
  EPS/SecMesh2D_fm.eps
  EPS/shape2.eps
  EPS/solve_build.eps
  EPS/THBaum.eps
  EPS/thinterpolate.eps
  EPS/turing_3D.eps
  EPS/turing.eps
  EPS/umalayout.eps)

set(SOURCES
  ../../src/course-examples/example01a_operator.hh
  ../../src/course-examples/example01a_Qk.hh
  ../../src/course-examples/example01a_RT.hh
  ../../src/course-examples/example01b_operator.hh
  ../../src/course-examples/example01b_Q2.hh
  ../../src/course-examples/example01.cc
  ../../src/course-examples/example01_main.hh
  ../../src/course-examples/example02_bcextension.hh
  ../../src/course-examples/example02_bctype.hh
  ../../src/course-examples/example02.cc
  ../../src/course-examples/example02_operator.hh
  ../../src/course-examples/example02_Q1.hh
  ../../src/course-examples/example03_bcextension.hh
  ../../src/course-examples/example03_bctype.hh
  ../../src/course-examples/example03.cc
  ../../src/course-examples/example03_operator.hh
  ../../src/course-examples/example03_Q2.hh
  ../../src/course-examples/example03_toperator.hh
  ../../src/course-examples/example04.cc
  ../../src/course-examples/example04.hh
  ../../src/course-examples/example04_operator.hh
  ../../src/course-examples/example05.cc
  ../../src/course-examples/example05_initial.hh
  ../../src/course-examples/example05_operator.hh
  ../../src/course-examples/example05_QkQk.hh
  ../../src/course-examples/example05_toperator.hh
  ../../src/course-examples/example06_bctype.hh
  ../../src/course-examples/example06.cc
  ../../src/course-examples/example06_Q1Q1.hh
  ../../src/course-examples/example07_adaptivity.hh
  ../../src/course-examples/example07.cc
  ../../src/course-examples/example07_error_indicator.hh
  ../../src/course-gridfunctionspace/analyticfunction.hh
  ../../src/course-gridfunctionspace/l2interpolationerror.hh
  ../../src/course-gridfunctionspace/parallelgfs.cc
  ../../src/course-gridfunctionspace/q1constrainedinterpolate.hh
  ../../src/course-gridfunctionspace/q1constrainedinterpolatemain.cc
  ../../src/course-gridfunctionspace/q1constraints.hh
  ../../src/course-gridfunctionspace/q1constraintsparameters.hh
  ../../src/course-gridfunctionspace/q1gridfunctionspace.hh
  ../../src/course-gridfunctionspace/q1gridfunctionspacemain.cc
  ../../src/course-gridfunctionspace/q1interpolate.hh
  ../../src/course-gridfunctionspace/q1interpolatemain.cc
  ../../src/course-gridfunctionspace/q1interpolationerror.hh
  ../../src/course-gridfunctionspace/q1interpolationerrormain.cc
  ../../src/course-gridfunctionspace/q1localbasis.hh
  ../../src/course-gridfunctionspace/q1localcoefficients.hh
  ../../src/course-gridfunctionspace/q1localfiniteelement.hh
  ../../src/course-gridfunctionspace/q1localfiniteelementmap.hh
  ../../src/course-gridfunctionspace/q1localinterpolation.hh
  ../../src/course-gridfunctionspace/q2interpolationerror.hh
  ../../src/course-gridfunctionspace/q2interpolationerrormain.cc
  ../../src/course-gridfunctionspace/thinterpolate.hh
  ../../src/course-gridfunctionspace/thinterpolatemain.cc
  ../../src/course-gridfunctionspace/thvelocity.hh
  src_examples/cadsample.cc
  src_examples/cadsample_operator.hh
  src_examples/cadsample_P1.hh
  src_examples/cadsample_parameter.hh
  src_examples/gridtest.cc

)

# Install documentation
# We assume that ${DOCFILES} and ${DOCFILES_EXTRA} are targets
if(LATEX_USABLE)
  dune_add_latex_document(pdelab-howto-article.tex FATHER_TARGET doc DEFAULT_PDF
      INPUTS ${TEX_FILES} turing.avi turing_3D.avi ${SOURCES} IMAGES ${EPS_IMAGES} BIBFILES lit.bib)
 create_doc_install(${DOCFILES} ${docdir} pdelab-howto-article_safepdf)
endif(LATEX_USABLE)
