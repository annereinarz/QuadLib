GCC=gcc
GCC_OPTS=-Wall -lm -g

all: 2dtria 2dquad

2dtria: squad2dtria

2dquad: squad2dquad

squad2dtria: SimpxSimp/squad2dtria.c  \
		SimpxSimp/squad2dtria_nomatlab.c  \
		SimpxSimp/sing_commonedge2dtria.c \
		SimpxSimp/sing_commonvertex2dtria.c \
		SimpxSimp/sing_distant2dtria.c \
		SimpxSimp/sing_identical2dtria.c \
		SimpxSimp/Ref2PhyS2d_affine.c \
		SimpxSimp/Ref2PhyS2d_spherical.c \
		quadpoints.c \
	   nomatlab_utils.c \
		mxMemory_nomatlab.c
	${GCC} ${GCC_OPTS} -o $@ $+

squad2dquad: CubexCube/squad2dquad.c  \
		CubexCube/squad2dquad_nomatlab.c  \
		CubexCube/sing_commonedge2dquad.c \
		CubexCube/sing_commonvertex2dquad.c \
		CubexCube/sing_identical2dquad.c \
		CubexCube/Ref2PhyJ2d.c \
		quadpoints.c \
        nomatlab_utils.c \
		mxMemory_nomatlab.c
	${GCC} ${GCC_OPTS} -o $@ $+
