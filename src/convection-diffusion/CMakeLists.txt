# $Id: duneproject 5425 2009-02-10 09:31:08Z sander $

set(SUBDIRS )

if(EIGEN_FOUND)
set(EIGENPROGRAMS
    laplacedirichletccfv-eigen)
endif(EIGEN_FOUND )

if(SUPERLU_FOUND)
set(SUPERLUPROGRAMS
    dgdiffusion
        dgdiffusion-hp
        dgldomain
        dgoverlapping
        rt0main)
endif(SUPERLU_FOUND )

set(noinst_PROGRAMS
    ${EIGENPROGRAMS}
    ${SUPERLUPROGRAMS}
        overlappingsinglephaseflowccfv
        overlappingsinglephaseflow
        nonoverlappingsinglephaseflow
        nonlineardiffusion
        transporttest
        instationarytest
        poisson
        laplacedirichletccfv
        reentrantcorner
        diffusion
        scalabilitytest
        ldomain)


set(noinst_HEADERS
        l2interpolationerror.hh
        parameterA.hh
        parameterB.hh
        parameterC.hh
        parameterD.hh
        parameterE.hh
        parameterF.hh
        problemA.hh
        problemB.hh
        problemC.hh
        problemD.hh
        problemE.hh
        problemF.hh)

set(dist_noinst_DATA
        grids/2dgrid.al
        grids/2dreentrantcorner.dgf
        grids/2dsimplex.alu
        grids/3dgrid.al
        grids/3drefinedcube.alu
        grids/CuKav.geo
        grids/CuKav.msh
        grids/circle.geo
        grids/circle1storder.msh
        grids/circle2ndorder.msh
        grids/cube.geo
        grids/cube.hexa
        grids/cube.tetra
        grids/cube1045.msh
        grids/cube3205.msh
        grids/cube99.msh
        grids/curved2d.geo
        grids/curved2d.msh
        grids/ldomain.al
        grids/ldomain.geo
        grids/ldomain.msh
        grids/ldomain2.msh
        grids/pyramid.geo
        grids/pyramid.msh
        grids/pyramid1storder.msh
        grids/pyramid2ndorder.msh
        grids/pyramid4.msh
        grids/telescope.geo
        grids/telescope.msh
        grids/telescope1storder.msh
        grids/telescope2ndorder.msh)

if(SUPERLU_FOUND)
add_executable(dgdiffusion dgdiffusion.cc)

add_dune_alberta_flags(dgdiffusion )
add_dune_amiramesh_flags(dgdiffusion )
add_dune_ug_flags(dgdiffusion )
add_dune_mpi_flags(dgdiffusion )
target_link_libraries(dgdiffusion ${DUNE_LIBS})
add_dune_superlu_flags(dgdiffusion )
add_dune_alugrid_flags(dgdiffusion )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgdiffusion_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgdiffusion_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(dgdiffusion ${DUNE_LIBS})
add_dune_superlu_flags(dgdiffusion )

add_executable(dgdiffusion_hp dgdiffusion-hp.cc)

add_dune_alberta_flags(dgdiffusion_hp )
add_dune_amiramesh_flags(dgdiffusion_hp )
add_dune_ug_flags(dgdiffusion_hp )
add_dune_mpi_flags(dgdiffusion_hp )
target_link_libraries(dgdiffusion_hp ${DUNE_LIBS})
add_dune_superlu_flags(dgdiffusion_hp )
add_dune_alugrid_flags(dgdiffusion_hp )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgdiffusion_hp_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgdiffusion_hp_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(dgdiffusion_hp ${DUNE_LIBS})
add_dune_superlu_flags(dgdiffusion_hp )

add_executable(dgldomain dgldomain.cc)
add_dune_alberta_flags(dgldomain )
add_dune_amiramesh_flags(dgldomain )
add_dune_ug_flags(dgldomain )
add_dune_mpi_flags(dgldomain )
target_link_libraries(dgldomain ${DUNE_LIBS})
add_dune_superlu_flags(dgldomain )
add_dune_alugrid_flags(dgldomain )
target_link_libraries(dgldomain ${DUNE_LIBS})
add_dune_superlu_flags(dgldomain )

add_executable(dgoverlapping dgoverlapping.cc)

add_dune_alberta_flags(dgoverlapping )
add_dune_amiramesh_flags(dgoverlapping )
add_dune_ug_flags(dgoverlapping )
add_dune_mpi_flags(dgoverlapping )
target_link_libraries(dgoverlapping ${DUNE_LIBS})
add_dune_superlu_flags(dgoverlapping )
add_dune_alugrid_flags(dgoverlapping )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgoverlapping_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(dgoverlapping_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})
endif(SUPERLU_FOUND )

add_executable(overlappingsinglephaseflowccfv overlappingsinglephaseflowccfv.cc)

add_dune_alberta_flags(overlappingsinglephaseflowccfv )
add_dune_amiramesh_flags(overlappingsinglephaseflowccfv )
add_dune_ug_flags(overlappingsinglephaseflowccfv )
add_dune_mpi_flags(overlappingsinglephaseflowccfv )
target_link_libraries(overlappingsinglephaseflowccfv ${DUNE_LIBS})
add_dune_superlu_flags(overlappingsinglephaseflowccfv )
add_dune_alugrid_flags(overlappingsinglephaseflowccfv )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(overlappingsinglephaseflowccfv_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(overlappingsinglephaseflowccfv_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(overlappingsinglephaseflowccfv ${DUNE_LIBS})
add_dune_superlu_flags(overlappingsinglephaseflowccfv )

add_executable(overlappingsinglephaseflow overlappingsinglephaseflow.cc)

add_dune_alberta_flags(overlappingsinglephaseflow )
add_dune_amiramesh_flags(overlappingsinglephaseflow )
add_dune_ug_flags(overlappingsinglephaseflow )
add_dune_mpi_flags(overlappingsinglephaseflow )
target_link_libraries(overlappingsinglephaseflow ${DUNE_LIBS})
add_dune_superlu_flags(overlappingsinglephaseflow )
add_dune_alugrid_flags(overlappingsinglephaseflow )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(overlappingsinglephaseflow_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(overlappingsinglephaseflow_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(overlappingsinglephaseflow ${DUNE_LIBS})
add_dune_superlu_flags(overlappingsinglephaseflow )

add_executable(nonoverlappingsinglephaseflow nonoverlappingsinglephaseflow.cc)

add_dune_alberta_flags(nonoverlappingsinglephaseflow )
add_dune_amiramesh_flags(nonoverlappingsinglephaseflow )
add_dune_ug_flags(nonoverlappingsinglephaseflow )
add_dune_mpi_flags(nonoverlappingsinglephaseflow )
target_link_libraries(nonoverlappingsinglephaseflow ${DUNE_LIBS})
add_dune_superlu_flags(nonoverlappingsinglephaseflow )
add_dune_alugrid_flags(nonoverlappingsinglephaseflow )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(nonoverlappingsinglephaseflow_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(nonoverlappingsinglephaseflow_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(nonoverlappingsinglephaseflow ${DUNE_LIBS})
add_dune_superlu_flags(nonoverlappingsinglephaseflow )

add_executable(poisson poisson.cc)

add_dune_alberta_flags(poisson )
add_dune_amiramesh_flags(poisson )
add_dune_ug_flags(poisson )
add_dune_mpi_flags(poisson )
target_link_libraries(poisson ${DUNE_LIBS})
add_dune_superlu_flags(poisson )
add_dune_alugrid_flags(poisson )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(poisson_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(poisson_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(poisson ${DUNE_LIBS})
add_dune_superlu_flags(poisson )

add_executable(nonlineardiffusion nonlineardiffusion.cc)

add_dune_alberta_flags(nonlineardiffusion )
add_dune_amiramesh_flags(nonlineardiffusion )
add_dune_ug_flags(nonlineardiffusion )
add_dune_mpi_flags(nonlineardiffusion )
target_link_libraries(nonlineardiffusion ${DUNE_LIBS})
add_dune_superlu_flags(nonlineardiffusion )
add_dune_alugrid_flags(nonlineardiffusion )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(nonlineardiffusion_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(nonlineardiffusion_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(nonlineardiffusion ${DUNE_LIBS})
add_dune_superlu_flags(nonlineardiffusion )

add_executable(transporttest transporttest.cc)

add_dune_alberta_flags(transporttest )
add_dune_amiramesh_flags(transporttest )
add_dune_ug_flags(transporttest )
add_dune_mpi_flags(transporttest )
message("transporttest DUNE_LIBS=${DUNE_LIBS}")
target_link_libraries(transporttest ${DUNE_LIBS})
add_dune_superlu_flags(transporttest )
add_dune_alugrid_flags(transporttest )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(transporttest_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(transporttest_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(transporttest ${DUNE_LIBS})
add_dune_superlu_flags(transporttest )

add_executable(instationarytest instationarytest.cc)

add_dune_alberta_flags(instationarytest )
add_dune_amiramesh_flags(instationarytest )
add_dune_ug_flags(instationarytest )
add_dune_mpi_flags(instationarytest )
target_link_libraries(instationarytest ${DUNE_LIBS})
add_dune_superlu_flags(instationarytest )
add_dune_alugrid_flags(instationarytest )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(instationarytest_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(instationarytest_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(instationarytest ${DUNE_LIBS})
add_dune_superlu_flags(instationarytest )

# mimetic is broken (probably irreparably...)

add_executable(laplacedirichletccfv laplacedirichletccfv.cc)

add_dune_ug_flags(laplacedirichletccfv )
add_dune_mpi_flags(laplacedirichletccfv )
target_link_libraries(laplacedirichletccfv ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(laplacedirichletccfv_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${UG_LDFLAGS} ${UG_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})

if(EIGEN_FOUND)
add_executable(laplacedirichletccfv_eigen laplacedirichletccfv.cc)


# The following lines are comment out as the translation is not clear. TODO review
#
# set(laplacedirichletccfv_eigen_ LDADD ${laplacedirichletccfv_LDADD})
endif(EIGEN_FOUND )

add_executable(reentrantcorner reentrantcorner.cc)

add_dune_alberta_flags(reentrantcorner )
add_dune_amiramesh_flags(reentrantcorner )
add_dune_ug_flags(reentrantcorner )
add_dune_mpi_flags(reentrantcorner )
add_dune_alugrid_flags(reentrantcorner )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(reentrantcorner_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(reentrantcorner_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(reentrantcorner ${DUNE_LIBS})
add_dune_superlu_flags(reentrantcorner )

if(SUPERLU_FOUND)
add_executable(rt0main rt0main.cc)

add_dune_alberta_flags(rt0main )
add_dune_amiramesh_flags(rt0main )
add_dune_ug_flags(rt0main )
add_dune_mpi_flags(rt0main )
target_link_libraries(rt0main ${DUNE_LIBS})
add_dune_superlu_flags(rt0main )
add_dune_alugrid_flags(rt0main )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(rt0main_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(rt0main_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})
endif(SUPERLU_FOUND )

add_executable(diffusion diffusion.cc)

add_dune_alberta_flags(diffusion )
add_dune_amiramesh_flags(diffusion )
add_dune_ug_flags(diffusion )
add_dune_mpi_flags(diffusion )
target_link_libraries(diffusion ${DUNE_LIBS})
add_dune_superlu_flags(diffusion )
add_dune_alugrid_flags(diffusion )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(diffusion_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(diffusion_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(diffusion ${DUNE_LIBS})
add_dune_superlu_flags(diffusion )

add_executable(scalabilitytest scalabilitytest.cc)

add_dune_alberta_flags(scalabilitytest )
add_dune_amiramesh_flags(scalabilitytest )
add_dune_ug_flags(scalabilitytest )
add_dune_mpi_flags(scalabilitytest )
target_link_libraries(scalabilitytest ${DUNE_LIBS})
add_dune_superlu_flags(scalabilitytest )
add_dune_alugrid_flags(scalabilitytest )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(scalabilitytest_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(scalabilitytest_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(scalabilitytest ${DUNE_LIBS})
add_dune_superlu_flags(scalabilitytest )

add_executable(ldomain ldomain.cc)

add_dune_alberta_flags(ldomain )
add_dune_amiramesh_flags(ldomain )
add_dune_ug_flags(ldomain )
add_dune_mpi_flags(ldomain )
target_link_libraries(ldomain ${DUNE_LIBS})
add_dune_superlu_flags(ldomain )
add_dune_alugrid_flags(ldomain )
# The following lines are comment out as the translation is not clear. TODO review
#
# set(ldomain_ LDADD ${DUNEMPILDFLAGS} ${DUNEMPILIBS} ${SUPERLU_LIBS} ${ALBERTA_LDFLAGS} ${ALBERTA_LIBS} ${AMIRAMESH_LDFLAGS} ${AMIRAMESH_LIBS} ${UG_LDFLAGS} ${UG_LIBS} ${ALUGRID_LDFLAGS} ${ALUGRID_LIBS} ${DUNE_LDFLAGS} ${DUNE_LIBS})
# The following lines are comment out as the translation is not clear. TODO review
#
# set(ldomain_ LDFLAGS ${AM_LDFLAGS} ${SUPERLU_LDFLAGS})

target_link_libraries(ldomain ${DUNE_LIBS})
add_dune_superlu_flags(ldomain )

set(CLEANFILES  poisson_*.vtu reentrant_*.vtu dgfparser.log q1*.vtu *.vtu *.png)

#============================================================
# pass most important options when "make distcheck" is used
set(DISTCHECK_CONFIGURE_FLAGS  --with-dune-common=${DUNE_COMMON_ROOT} --with-dune-grid=${DUNE_GRID_ROOT} --with-dune-istl=${DUNE_ISTL_ROOT} --with-dune-localfunctions=${DUNE_LOCALFUNCTIONS_ROOT} --with-dune-pdelab=${DUNE_PDELAB_ROOT}  CXX="${CXX}" CC="${CC}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules

foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
