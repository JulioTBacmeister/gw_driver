##################################
#   $
##################################

VPATH=gwdrag/gw_drag

XESM_SRCS = \
	shr_kind_mod.F90   \
	shr_const_mod.F90   \
	shr_infnan_mod.F90   \
	shr_strconvert_mod.F90   \
	shr_log_mod.F90   \
	shr_sys_mod.F90   \
	dycore_mod.F90 \
	coords_1d.F90 \
	gw_utils.F90        \
	interpolate_data.F90 \
	linear_1d_operators.F90 \
	local_physconst.F90        \
	local_constituent.F90        \
	vdiff_lu_solver.F90        \
	gw_diffusion.F90        \
	gw_common.F90        \
	local_physics_types.F90        \
	gw_convect.F90        \
	gw_rdg.F90        \
	xesm.F90


COMPILER=gfortran
XESM_OBJS= ${XESM_SRCS:.F90=.o}
DRYTEST_OBJS= ${DRYTEST_SRCS:.F90=.o}
ONEDTEST_OBJS= ${ONEDTEST_SRCS:.F90=.o}

OSYST := $(shell uname)



#$(OBJS) : %.o: %.F90
#XXgoldyXX: Changed this into an implicit rule to allow VPATH to work
%.o: %.F90
	$(COMPILER) $(F90_FLAGS) $<
#	$(COMPILER) $(F90_FLAGS) $*.F90


#STD_FLAGS=-c
STD_FLAGS=-c -g -fbacktrace -fdollar-ok -fcheck='all' -ffpe-trap='invalid','zero','overflow','underflow' -ffpe-summary='all' -ffree-line-length-0

# -trap=INVALID

machine=my_mac
#ifeq ($(machine),my_mac)
    INC_NETCDF := /Users/juliob/opt/local/include
    LIB_NETCDF := /Users/juliob/opt/local/lib
#  else

#    INC_NETCDF :=/usr/local/netcdf-c-4.6.1-f-4.4.4-gcc-g++-gfortran-4.8.5/include
#    LIB_NETCDF :=/usr/local/netcdf-c-4.6.1-f-4.4.4-gcc-g++-gfortran-4.8.5//lib
#  endif

    FFLAGS   := -I$(INC_NETCDF)


  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lnetcdff
  #FFLAGS   := -c -Mlarge_arrays -I$(INC_NETCDF)
  #ifeq ($(OSYST),Linux)
  #  FFLAGS   := -I$(INC_NETCDF)
  #endif
  #ifeq ($(OSYST),Darwin)
  #  FFLAGS   := -I$(INC_NETCDF)
  #endif



F90_FLAGS=$(FFLAGS) $(STD_FLAGS) $(USER_FDEFS) -DUSE_CONTIGUOUS= -DCPRGNU -DUNITTEST


CPP_FLAGS=-P -DALPHA_MACH



#-----------------------------------------------------
# These files (RHS) are used by almost everyone.
# So, safer just to rebuild everything if they change.
#-----------------------------------------------------

report:
	echo $(machine)

uxesm  : $(XESM_OBJS)
	$(COMPILER) -O1 $(XESM_OBJS) $(LDFLAGS) -o uxesm.x

clean : 
	/bin/rm -f *.o *.a *.x *.mod fort.511

reset : 
	/bin/rm -f *.o *.a *.x *.mod gw_rdg_calc_inc.F90 \
	gw_drag.F90 gw_common.F90 gw_diffusion.F90 gw_rdg.F90 \
	gw_utils.F90 gw_oro.F90 vdiff_lu_solver.F90 \
	linear_1d_operators.F90 coords_1d.F90 physconst.F90 \
	shr_kind_mod.F90 shr_const_mod.F90 shr_flux_mod.F90 \
	shr_log_mod.F90 shr_strconvert_mod.F90
