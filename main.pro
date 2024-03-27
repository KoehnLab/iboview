# Default compilation command sequence:
#
#     mkdir build && cd build
#     qmake-qt4 ../main.pro && make -j 12 && make install
#     cd ..
#
# Afterwards there should be a "iboview" executable in the main directory ("make
# install" copies the exe from the build directory to the project file
# directory)
#
# Note: if the external blas does not work, you can compile with
#
#     BLASLIB= MKLROOT= qmake-qt4 main.pro && make -j 12 && make install
#
# to force internal LAPACK emulation

############################### CONFIG: Hacky overrides for cgk's computer

win32 {
#   BOOST_ROOT = C:\SDKs\boost_1_77_0   # <-- yes, that is where it is on my own computer.
}

############################### CONFIG: Check for availability of MKL or other
############################### external BLAS/LAPACK

BLASLIB=$$(BLASLIB)
# ^-- pick up environment variable BLASLIB if present. See: https://doc.qt.io/archives/qt-4.8/qmake-advanced-usage.html
!isEmpty(BLASLIB) {
   message(using BLAS linker line from environment variable BLASLIB: $$BLASLIB)
}
# no BLAS/LAPACK linker line explicitly specified as environment variable?
isEmpty(BLASLIB) {
   # try to pick up MKL from the (default) MKLROOT environment variable.
   # You may need to "load module mkl" before doing that.
   # If the MKLROOT environment variable is not set, you can also override
   # the path manually (e.g., setting MKLDIR = /opt/intel/composerxe/mkl)
   MKLDIR = $$(MKLROOT)
   !isEmpty(MKLDIR) {
      exists($$MKLDIR) {
         message(Using MKL from MKLROOT: $$MKLDIR)
         BLASLIB = -L$$MKLDIR/lib/intel64 -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -fopenmp -lpthread -Wl,-rpath,$$MKLDIR/lib/intel64
      }
      !exists($$MKLDIR) {
         message(Directory $$MKLDIR specified in MKLROOT environment variable does not exist)
      }
   }
   # To set up some other sort of BLAS/LAPACK manually: set the required linker
   # line to BLASLIB. E.g.:
   #
   # BLASLIB = -LC:\Users\cgk\Develop\Libraries\OpenBLAS\bin -lopenblas
   #
   # Notes: 
   # - Header or include files are not needed, and will not be used if
   #   provided.
   # - Make sure that is a BLAS with 64bit integer interface!
   #   (e.g., -lmkl_intel_ilp64, NOT -lmkl_intel_lp64).
}

!isEmpty(BLASLIB) {
   message(Using external BLAS/LAPACK: $$BLASLIB)
   LIBS += $$BLASLIB
   HAVE_EXTERNAL_BLAS = 1
}
isEmpty(BLASLIB) {
   message(Using internal BLAS/LAPACK emulation.)
   # set HAVE_EXTERNAL_BLAS to 0. In this case, we will include compilation of
   # our hacked f2c LAPACK and a eigen BLAS wrapper. This is normally good
   # enough for IboView analysis applications (MicroScf might suffer a bit,
   # though)
   HAVE_EXTERNAL_BLAS = 0
}

############################### CONFIG: Check for C++ boost libraries
############################### (we only use header-only sub-libraries)

isEmpty(BOOST_ROOT) {          # <-- boost directory not manually set in .pro file?
   BOOST_ROOT=$$(BOOST_ROOT)   # <-- check if we got it from an environment variable
}
!isEmpty(BOOST_ROOT) : exists($$BOOST_ROOT/boost/math/special_functions/beta.hpp) {
   CONFIG += config_boost_special  # found it. Don't do further tests.
}
!lessThan(QT_MAJOR_VERSION, 5) : !CONFIG(config_boost_special) {
   message(boost test: qtCompileTest branch)
   load(configure) # to enable qtCompileTest (unfortunately, that's QT5+ only...)

   # Enable test for presence of boost libraries (optional).
   # - https://stackoverflow.com/questions/22238446/qmake-checking-for-library
   # - https://doc.qt.io/qt-5/qmake-test-function-reference.html#qtcompiletest-test
   qtCompileTest(boost_special)
}
lessThan(QT_MAJOR_VERSION, 5) : !CONFIG(config_boost_special) {
   # try to pick it up at a default location
   isEmpty(BOOST_ROOT) {
      BOOST_ROOT=/usr/include
   }
}
!isEmpty(BOOST_ROOT) : exists($$BOOST_ROOT/boost/math/special_functions/beta.hpp) {
   CONFIG += config_boost_special
   INCLUDEPATH += $$BOOST_ROOT
   message(Found boost::math special functions at $$BOOST_ROOT/boost/math/special_functions)
   # enabling this will enable the log Gamma optional functionals for IBO localization
   # (which use gamma/digamma/trigamma functions) and a cgk special for radial DFT grids
   # using incomplete beta functions. Boost special functions are header only, and do not
   # need compiled libraries. However, they are still a pain to install on some platforms,
   # so this functionality is optional.
   DEFINES += HAVE_BOOST_SPECIAL_FUNCTIONS
}
#!CONFIG(config_boost_special) {
#   message(Disabling optional features based on boost/math/special_functions (not found))
#}


###! START_DIRS
MICROSCF = src/MicroScf
MIGRID = src/MicroScf
IRCORE = src/IrCore
CXF = src/Common
IV = src/IboView
OTHER = src/External
GL = src/GL
PUGIXML = src/pugixml
equals(HAVE_EXTERNAL_BLAS,0) {
   EIGEN = src/eigen
   EIGENBLAS = src/eigen/blas
   INCLUDEPATH += $$EIGEN
   INCLUDEPATH += $$EIGENBLAS
}
###! END_DIRS

TEMPLATE = app
TARGET = iboview
DEPENDPATH += .
INCLUDEPATH += src
INCLUDEPATH += $$IV
INCLUDEPATH += $$MICROSCF
INCLUDEPATH += $$IRCORE
INCLUDEPATH += $$CXF
INCLUDEPATH += $$OTHER
INCLUDEPATH += $$PUGIXML
DEFINES += USE_CTINT1E_H
DEFINES += PROG_IBOVIEW
DEFINES += ENABLE_GLOBAL_RNG
DEFINES += IR_ECP
DEFINES += IR_KERNEL_PTRS
DEFINES += CX_ASSERT


# By default, make "make install" copy the final executable to the same
# directory in which the main project file (main.pro) lies. IboView has all its
# data files embedded as QT resources, so apart from the exe itself no extra
# files are needed.
target.path = $$_PRO_FILE_PWD_
# program.files = $$TARGET
INSTALLS += target
# ^-- see: https://doc.qt.io/qt-5/qmake-advanced-usage.html#installing-files


# LIBS += -lGLEW     # OpenGL extension wrangler library (rest is included automatically via QT += opengl).
# (brings own implementation of glew.c now, from glew-20140726.tgz)
DEFINES += GLEW_STATIC
SOURCES += $$GL/glew.c
HEADERS += $$GL/glew.h $$GL/glxew.h $$GL/wglew.h

QT += core gui widgets opengl script svg
CONFIG += qt
CONFIG -= debug    # turn off debug build
CONFIG += release  # turn on release build
DEFINES += NDEBUG  # turn of assertions

# OpenMP
win32 {
   QMAKE_CXXFLAGS += /openmp
   QMAKE_LFLAGS += /openmp
}
!win32 {
   QMAKE_CXXFLAGS += -fopenmp
   QMAKE_LFLAGS += -fopenmp
}
# ^- warning: need to replace all OpenMP loop variables by "int". Otherwise won't work
#    in vc or gcc <= 4.4 (unsigned integral types only allowed since OpenMP 3.0... wtf
#    were they thinking?!)

!win32 {
   # Flags for g++ or clang (in g++ emulation mode on linux)
   OPTFLAGS_OFF = -O2
   OPTFLAGS_ON = -Ofast -ffast-math -march=native
   QMAKE_CXXFLAGS += -Wno-deprecated-copy -Wno-class-memaccess -Wno-implicit-fallthrough -Wno-unknown-pragmas -fmax-errors=10
   # DEFINES += _GLIBCXX_PARALLEL
   # ^-- Interesting, see: https://gcc.gnu.org/onlinedocs/libstdc++/manual/parallel_mode_using.html
   #     But still not sure what to think about this. In any case, we
   #     do not use any of those algorithms in places where it would help.
}
win32 {
   # Flags for MSVC on Windows (Note: "win32" also includes 64bit configurations)
   # (we assume we use MSVC on Windows. The actual compiler specialization
   # mechanism for some reason does not work on my qmake).
   QMAKE_CXXFLAGS += /MP
   LIBS += opengl32.lib
   RC_FILE = resources/windows_icon_v2.rc
   OPTFLAGS_OFF =
   OPTFLAGS_ON = /MP /O2 /openmp /arch:AVX
   TARGET = IboView
   target.path = $$_PRO_FILE_PWD_/bin
}
QMAKE_CFLAGS_RELEASE  -= $$OPTFLAGS_OFF
QMAKE_CFLAGS_RELEASE  += $$OPTFLAGS_ON
QMAKE_CXXFLAGS_RELEASE -= $$OPTFLAGS_OFF
QMAKE_CXXFLAGS_RELEASE += $$OPTFLAGS_ON

# Include experimental/in-development/testing code. Note that most such code is
# actually filtered out by my deploy scripts in release versions of the code,
# so you may not actually see any of that in a non-upstream version of the
# programs
#DEFINES += INCLUDE_OPTIONALS


# Input
FORMS += $$IV/MainForm2.ui $$IV/AboutForm.ui $$IV/FindOrbitalsForm.ui $$IV/ShowTextForm.ui $$IV/ComputeWfForm.ui $$IV/TablesForm.ui $$IV/EditFramesForm.ui $$IV/PreferencesForm.ui $$IV/ComputeEosForm.ui $$IV/$$IV/EditVolumeSurfaceForm.ui
RESOURCES += resources.qrc

equals(HAVE_EXTERNAL_BLAS,0) {
   # No BLAS/LAPACK set up :(. I here follow the philosophy that that the most
   # important performance improvement a program can undergo is the transition
   # from the non-working into the working state. So I spent quite some time on
   # creating workarounds to arrive at a program which may not run 100% optimal,
   # but will work:
   #
   # - I reduced the amount of different LAPACK-routines used rather significantly
   #   throughout both MicroScf and IboView. And the few routines which are
   #   remaining have been programmed in such a way that they can perform their
   #   tasks with LAPACK routines replaced by Eigen routines.
   #   At this moment, the only higher level LA primitives used are:
   #
   #   + Diagonalize(): Square symmetric matrix spectral decomposition
   #     (H -> H = V diag(lambda) V^T):
   #     Normally implemented with LAPACK DSYEVD (divide&conquor) or DSYEVR (RRR).
   #     In USE_EIGEN_LA mode, this is replaced by a implementation using
   #     Eigen::SelfAdjointEigenSolver (Householder QR) on Eigen::Map<>s. That
   #     has some performance issues, especially in terms of parallel scaling (see below),
   #     but should be rather robust and accurate.
   #
   #   + ComputeSvd(): General singular value decomposition (A -> A = U diag(sigma) V.T):
   #     Normally implemented with LAPACK DGESDD (divide&conquor). In USE_EIGEN_LA mode,
   #     replaced with a implementation based on Eigen::JacobiSVD<> on Eigen::Map<>s.
   #     For SVDs, robustness and accuracy are the core criteria, and we do not use them
   #     anywhere where performance is important. I think the Eigen implementation for
   #     those should be a suitable replacement in allmost all our applications.
   #
   #   + CalcCholeskyFactors(), TriangularSolve(), TriangularMxm(), TriangularInvert():
   #     These are really more BLAS-like than LA-like (very simple and straight-forward
   #     formal implementations, and the core question is how to optimize it for a platform).
   #     They are replaced by Eigen:LLT on Eigen::Map<>s.
   #
   # - We will replace the BLAS with Eigen's BLAS emulation.
   #
   #   + It comes with Eigen's own distribution since quite a while. I just copied
   #     and linked in the Eigen release as-is
   #
   #   + The Eigen-BLAS is not actually that bad: the dense algebra kernels
   #     which are well implemented can make efficient use of the CPU thread they
   #     are running on, with a rather nice AVX/AVX2 (in the meantime) kernel setup.
   #     I've seen comparable performance to MKL in many situations.
   #
   #   + ...However, that is the case provided it is compiled with a modern optimizer
   #     for numerical work (which we will not) and for the concrete CPU at han
   #     (which we also will not do). We just make a static AVX version -- so
   #     all routines which rely cruvially on efficient MxMs (e.g., DF-JK) will
   #     be several times slower than what they could achieve with AVX2, FMA,
   #     AVX512, etc. For DF-RKS (especially in DF-JX mode) the effect is not so
   #     pronounced, however, and it can run quite acceptably if compiled like we do here
   #
   #   + The most serious issue with Eigen-BLAS (and the emulated LAPACK) is the
   #     nearly complete lack of multi-threading parallelization. At the time of writing,
   #     the only Eigen kernel with any degree of OpenMP parallelization is the
   #     dense-square-matrix MxM core, and even that one is not done particularly well.
   #     Several other routines which are crucial for many QC algorithms (e.g.,
   #     triangular mxm/solve, Cholesky decomposition, QR decomposition&co) are quite
   #     efficient if run on a single core, but have no parallelization whatsoever.
   #     And some routines which are very helpful in special cases are just missing
   #     (e.g., symmetric rank-K updates (syrk) --- effectively matrix multiplications
   #     A * A.T which are known to yield a symmetric product. Saves factor 2 over
   #     straight MxM).
   #
   #   + For most tasks IboView would be used for those things are not really a
   #     deal-breaker, but even here such LA primitives are sprinkled all throughout
   #     the program, and those not being paralellized at all puts a rather tight
   #     upper limit on what can be achieved in terms of parallel scaling. E.g., this
   #     leads to the Fock matrix diagonalization have zero parallel scaling in the
   #     integrated SCF routines, and on faster machines and larger systems this will
   #     therefore quickly dominate the total cost of DF-RKS.
   #
   SOURCES += $$EIGENBLAS/double.cpp $$EIGENBLAS/xerbla.cpp
   DEFINES += NO_OVERWRITE USE_EIGEN_LA
}

HEADERS += $$IV/Iv.h $$IV/IvMain.h $$IV/IvSettings.h $$IV/IvView3D.h $$IV/IvDataSet.h $$IV/IvDataOptions.h $$IV/IvVolumeDataSet.h $$IV/IvOrbital.h $$IV/IvDocument.h $$IV/IvAnalysis.h $$IV/IvIsoSurface.h $$IV/IvMesh.h $$IV/IvScript.h $$IV/IvIao.h $$IV/IvIrc.h $$IV/IvGl.h $$IV/IvLog.h $$IV/fn_LiberationSans.h
SOURCES += $$IV/IvMain.cpp $$IV/IvSettings.cpp $$IV/IvView3D.cpp $$IV/IvDataSet.cpp $$IV/IvDataOptions.cpp $$IV/IvVolumeDataSet.cpp $$IV/IvOrbital.cpp $$IV/IvDocument.cpp  $$IV/IvAnalysis.cpp $$IV/IvScript.cpp $$IV/IvMesh.cpp $$IV/IvIsoSurface.cpp $$IV/IvIao.cpp $$IV/IvIrc.cpp $$IV/IvGl.cpp $$IV/IvLog.cpp

HEADERS += $$IV/IvOrbitalFile.h $$IV/IvCurveView.h $$IV/IvComputeWfForm.h $$IV/IvFindOrbitalsForm.h $$IV/IvShowTextForm.h $$IV/IvFixedAspectSvg.h $$IV/IvStatusBar.h $$IV/IvTables.h $$IV/IvEditFramesForm.h $$IV/IvPreferencesForm.h $$IV/IvEditVolumeSurface.h $$IV/IvFileConvert.h $$IV/IvComputeEosForm.h
SOURCES += $$IV/IvOrbitalFile.cpp $$IV/IvCurveView.cpp $$IV/IvComputeWfForm.cpp $$IV/IvFindOrbitalsForm.cpp $$IV/IvShowTextForm.cpp $$IV/IvFixedAspectSvg.cpp $$IV/IvStatusBar.cpp $$IV/IvTables.cpp $$IV/IvEditFramesForm.cpp $$IV/IvPreferencesForm.cpp $$IV/IvEditVolumeSurface.cpp $$IV/IvFileConvert.cpp $$IV/IvComputeEosForm.cpp
# SOURCES += $$IV/IvIntInit.cpp 


HEADERS += $$PUGIXML/pugixml.hpp $$PUGIXML/pugiconfig.hpp
SOURCES += $$PUGIXML/pugixml.cpp

SOURCES += $$OTHER/memory_size.c
HEADERS += $$OTHER/optionparser.h
HEADERS += $$OTHER/vector_math.h

HEADERS += $$OTHER/QPropertyModel.h
SOURCES += $$OTHER/QPropertyModel.cpp


HEADERS += $$IRCORE/Ir.h
HEADERS += $$IRCORE/IrFactorials.inl
HEADERS += $$IRCORE/IrFactorials.h
HEADERS += $$IRCORE/IrAmrr.h
SOURCES += $$IRCORE/IrAmrr.cpp
HEADERS += $$IRCORE/IrBoysFn.h
SOURCES += $$IRCORE/IrBoysFn.cpp
HEADERS += $$IRCORE/IrBoysFn.inl
SOURCES += $$IRCORE/IrCore.cpp
SOURCES += $$IRCORE/IrDrv.cpp
SOURCES += $$IRCORE/IrGridOps.cpp
HEADERS += $$IRCORE/IrInternal.h
HEADERS += $$IRCORE/IrDrvAux.h
SOURCES += $$IRCORE/IrDrvAux.cpp
SOURCES += $$IRCORE/IrMeta2i.cpp
HEADERS += $$IRCORE/IrMeta2i.h
SOURCES += $$IRCORE/IrSlmRot.cpp
SOURCES += $$IRCORE/IrImportTrafo.cpp
HEADERS += $$IRCORE/IrEcp.h
SOURCES += $$IRCORE/IrEcp.cpp
SOURCES += $$IRCORE/IrEcpBesselFn.cpp
SOURCES += $$IRCORE/IrComplexRlm.cpp
HEADERS += $$IRCORE/IrComplexRlm.h
HEADERS += $$IRCORE/IrComplexRlm.h
# SOURCES += $$IRCORE/IrSlmX.cpp
# ^- part of IrAmrr.cpp now


HEADERS += $$CXF/CxAngularGrid_Orbits.inl
HEADERS += $$CXF/CxAngularGrid_Grids.inl
HEADERS += $$CXF/CxAngularGrid.h
SOURCES += $$CXF/CxAngularGrid.cpp
DEFINES += USE_AIGG_GRIDS
# SOURCES += $$IV/CtDftGrid_ivb.cpp
# HEADERS += $$IV/CtDftGrid_ivb.h
HEADERS += $$MIGRID/CtDftGrid.h
SOURCES += $$MIGRID/CtDftGrid.cpp
HEADERS += $$MIGRID/CtDftGrid_Params.h
SOURCES += $$MIGRID/CtDftGrid_Params.cpp
HEADERS += $$MIGRID/CtDftGrid_Radial.h
SOURCES += $$MIGRID/CtDftGrid_Radial.cpp
SOURCES += $$MIGRID/CtDftGrid_CommonDefs.h
HEADERS += $$MIGRID/CtDftGrid_QuadCriteria.h
SOURCES += $$MIGRID/CtDftGrid_QuadCriteria.cpp
SOURCES += $$MIGRID/CtVoronoiPartition.cpp
HEADERS += $$MIGRID/CtVoronoiPartition.h
SOURCES += $$CXF/CxDensityModel.cpp
HEADERS += $$CXF/CxDensityModel.h
HEADERS += $$MICROSCF/CtAtomDensity.inl
HEADERS += $$MICROSCF/CtAtomDensity.h
SOURCES += $$MICROSCF/CtAtomDensity.cpp
# HEADERS += atomic_density_tab_backup/CtAtomDensity.inl
# HEADERS += atomic_density_tab_backup/CtAtomDensity.h
# SOURCES += atomic_density_tab_backup/CtAtomDensity.cpp

HEADERS += $$MICROSCF/CtAtomSet.h
SOURCES += $$MICROSCF/CtAtomSet.cpp
HEADERS += $$MICROSCF/CtBasisLibrary.h
SOURCES += $$MICROSCF/CtBasisLibrary.cpp
HEADERS += $$MICROSCF/CtBasisSet.h
SOURCES += $$MICROSCF/CtBasisSet.cpp
HEADERS += $$MICROSCF/CtBasisShell.h
SOURCES += $$MICROSCF/CtBasisShell.cpp
HEADERS += $$MICROSCF/CtBasisDesc.h
SOURCES += $$MICROSCF/CtBasisDesc.cpp
HEADERS += $$MICROSCF/CtInt1e.h
SOURCES += $$MICROSCF/CtInt1e.cpp
# SOURCES += $$MICROSCF/CtMain.cpp
SOURCES += $$MICROSCF/CtMatrix.cpp
HEADERS += $$MICROSCF/CtMatrix.h
HEADERS += $$MICROSCF/CtWfi.h
SOURCES += $$MICROSCF/CtWfi.cpp
HEADERS += $$MICROSCF/CtWfa.h
SOURCES += $$MICROSCF/CtWfa.cpp
HEADERS += $$MICROSCF/CtIvStubs.h
# ^-- note: the .cpp of it should *not* be included!
HEADERS += $$MICROSCF/CtVolumeProperty.h
SOURCES += $$MICROSCF/CtVolumeProperty.cpp
SOURCES += $$MICROSCF/CtAtomConfigs.cpp
SOURCES += $$MICROSCF/CtDft.cpp
HEADERS += $$MICROSCF/CtDft.h
SOURCES += $$MICROSCF/CtDfti.cpp
HEADERS += $$MICROSCF/CtDfti.h
SOURCES += $$MICROSCF/CtDftFunc.cpp
HEADERS += $$MICROSCF/CtDftFunc.h
SOURCES += $$MICROSCF/CtFockBuild.cpp
HEADERS += $$MICROSCF/CtFockBuild.h
SOURCES += $$MICROSCF/CtRhf.cpp
HEADERS += $$MICROSCF/CtRhf.h
HEADERS += $$MICROSCF/CtRhfOptions.h
SOURCES += $$MICROSCF/CtRhfOptions.cpp
SOURCES += $$MICROSCF/CtRhfGuess.cpp
HEADERS += $$MICROSCF/CtRhfGuess.h
SOURCES += $$MICROSCF/CtIao.cpp
HEADERS += $$MICROSCF/CtIao.h
SOURCES += $$MICROSCF/CtOrbLoc.cpp
HEADERS += $$MICROSCF/CtOrbLoc.h
HEADERS += $$MICROSCF/CtWfDecl.h
SOURCES += $$MICROSCF/CtWfDecl.cpp
SOURCES += $$MICROSCF/xc_meta.cpp
HEADERS += $$MICROSCF/xc_meta.h
SOURCES += $$MICROSCF/coords_meta.cpp
HEADERS += $$MICROSCF/coords_meta.h
HEADERS += $$MICROSCF/CtDftDispCorr.h
HEADERS += $$MICROSCF/CtDftDispCorr.inl
SOURCES += $$MICROSCF/CtDftDispCorr.cpp
# HEADERS += $$MICROSCF/CtCosmo.h
# SOURCES += $$MICROSCF/CtCosmo.cpp

HEADERS += $$MICROSCF/format.h
SOURCES += $$MICROSCF/format.cpp
HEADERS += $$CXF/CxArgSort.h
HEADERS += $$CXF/CxParse1.h
SOURCES += $$CXF/CxParse1.cpp
SOURCES += $$CXF/CxAssertFail.cpp
SOURCES += $$CXF/CxIntrusivePtr.h
HEADERS += $$CXF/CxRawAtom.h
SOURCES += $$CXF/CxRawAtom.cpp
HEADERS += $$CXF/CxXyzFrameIo.h
SOURCES += $$CXF/CxXyzFrameIo.cpp
HEADERS += $$CXF/CxIo.h
SOURCES += $$CXF/CxIo.cpp
HEADERS += $$CXF/CxTiming.h
SOURCES += $$CXF/CxTiming.cpp
HEADERS += $$CXF/CxPhysicalUnits.h
HEADERS += $$CXF/CxAtomData.h
HEADERS += $$CXF/CxAtomData.inl
SOURCES += $$CXF/CxAtomData.cpp
HEADERS += $$CXF/CxSortedTable.h
HEADERS += $$CXF/CxSequenceIo.h
HEADERS += $$CXF/CxIterTools.h
HEADERS += $$CXF/CxPrintLevel.h
SOURCES += $$CXF/CxPrintLevel.cpp
HEADERS += $$CXF/CxAtomParamSpec.h
HEADERS += $$CXF/CxNewtonAux.h
HEADERS += $$CXF/CxNewtonAux.cpp
HEADERS += $$CXF/CxAlgebra.h
SOURCES += $$CXF/CxAlgebra.cpp
HEADERS += $$CXF/CxDefs.h
HEADERS += $$CXF/CxDiis.h
SOURCES += $$CXF/CxDiis.cpp
HEADERS += $$CXF/CxFortranInt.h
HEADERS += $$CXF/CxMemoryStack.h
SOURCES += $$CXF/CxMemoryStack.cpp
HEADERS += $$CXF/CxOpenMpProxy.h
HEADERS += $$CXF/CxOpenMpAcc.h
HEADERS += $$CXF/CxPodArray.h
HEADERS += $$CXF/CxStorageDevice.h
SOURCES += $$CXF/CxStorageDevice.cpp
HEADERS += $$CXF/CxTypes.h
HEADERS += $$CXF/CxVec3.h
# HEADERS += $$CXF/CxMath.h
HEADERS += $$CXF/CxMathAlgo.h
HEADERS += $$CXF/CxColor.h
SOURCES += $$CXF/CxColor.cpp
HEADERS += $$CXF/CxPoly.h
SOURCES += $$CXF/CxPoly.cpp
HEADERS += $$CXF/CxOsInt.h
SOURCES += $$CXF/CxOsInt.cpp
# HEADERS += $$CXF/CxMiniPointGroup.h
# SOURCES += $$CXF/CxMiniPointGroup.cpp
# HEADERS += $$CXF/CxBitOps.h
# SOURCES += $$CXF/CxBitOps.cpp
HEADERS += $$CXF/CxRandom.h
SOURCES += $$CXF/CxRandom.cpp
# ^--CxRandom used for random rotation code in orbital localization to test for stability
#    of finding global optimum.
# HEADERS += $$CXF/CxIndentStream.h
# SOURCES += $$CXF/CxIndentStream.cpp
# SOURCES += $$CXF/CxNumpyArray.cpp
# HEADERS += $$CXF/CxNumpyArray.h



# https://stackoverflow.com/questions/18371797/qmake-handling-options-for-both-gcc-and-msvc
# gcc:QMAKE_CXXFLAGS += -Wno-deprecated-copy -Wno-class-memaccess
# message($${QMAKE_CXXFLAGS})
# message($${QMAKE_COMPILER})
# ^- doesn't work. Apparently no compiler is set.


#CONFIG -= dynamic
#CONFIG += static
# ^- I wonder if that works? ... nop. no effect at all. Still need 10324324432 QT DLLs,
#    and in the very same directory as the exe file.

# note: for compiling with clang instead of g++ on linux:
#   /usr/bin/qmake-qt4 -spec /usr/lib64/qt4/mkspecs/unsupported/linux-clang
# (works for me, but doesn't do OpenMP)



# Config testing:
#
#    rm -f .qmake.cache .qmake.stash && qmake-qt5 main.pro
#    rm -f .qmake.cache .qmake.stash && qmake-qt4 main.pro
#    rm -f .qmake.cache .qmake.stash && BOOST_ROOT=/usr/include qmake-qt4 main.pro
#    rm -f .qmake.cache .qmake.stash && BOOST_ROOT=/certainly-not-here qmake-qt4 main.pro
#
# Expected results:
# - First one should run config.test/boost_special and, depending on what it finds there, pass (works for me)
# - Second should pass (find boost in default location at /usr/include/boost)
# - Third one should pass (directory explicitly specified, and found the header there)
# - Fourth one should fail (no boost library to be found there)
#
# Note: need to delete cache files to make qmake retest the configuration.
# Otherwise it sometimes comes to rather odd and incoherent results if configs
# are mixed.
