GCC=gcc
GCC_OPTS=-Wall -lm -g

all: 1d

1d: squad1d

squad1d:                squad1d_nomatlab.c           \
                            squad1d.c                   \
                            sing_identical1d.c       \
                            sing_commonvertex1d.c    \
                            Ref2PhyS1d_affine.c      \
                            Ref2PhyS1d_circle.c      \
                            quadpoints.c \
									 nomatlab_utils.c \
									 mxMemory_nomatlab.c
	${GCC} ${GCC_OPTS} -o $@ $+
