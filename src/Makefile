##
## SWD-ELL trans-dimensional inversion makefile (RF/MT machinery removed)
## Jan Dettmer jand@uvic.ca 2012; SWD-ELL-PMPI 2026
##
include ./Makefile.compiler
##############################################################################
BIN = ./bin
OBJ = ./obj
OBJS = $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/alloc_obj.o $(OBJ)/read_input.o $(OBJ)/loglhood.o \
       $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/quicksort.o \
       $(OBJ)/UPDATE_MAPfile.o $(OBJ)/xcorr.o $(OBJ)/COVmatEst.o $(OBJ)/UpdateCOV.o \
       $(OBJ)/prjmh_temper_rf.o
FO = $(OBJ)/libswd.a
LIB = -framework Accelerate   # LAPACK/BLAS for DPOTRF/DPOTRI (COVmatEst); use -llapack -lblas on Linux
##############################################################################

all: dirs bin/prjmh_temper_rf

dirs:
	@mkdir -p $(BIN) $(OBJ)

prjmh_temper_rf: bin/prjmh_temper_rf

clean:
	rm -f bin/* $(OBJ)/* *.o *.mod
	(cd ./swd/;/bin/rm -f *.o *.mod *.l core *.a;)

$(BIN)/prjmh_temper_rf: libswd $(OBJS) $(FO)
	$(FC)  -o $(BIN)/prjmh_temper_rf $(OBJS) $(FO) $(LIB)

$(OBJ)/prjmh_temper_rf.o: prjmh_temper_rf.f90 $(OBJ)/rjmcmc_com.o $(OBJ)/nr.o
	$(FC) -c -o $(OBJ)/prjmh_temper_rf.o prjmh_temper_rf.f90

$(OBJ)/rjmcmc_com.o: rjmcmc_com.f90 $(OBJ)/data_type.o
	$(FC) -c rjmcmc_com.f90 -o $(OBJ)/rjmcmc_com.o

$(OBJ)/data_type.o: data_type.f90
	$(FCSERIAL) -c data_type.f90 -o $(OBJ)/data_type.o

libswd:
	(cd ./swd/;make;)
	/bin/ln -sf ./swd/*.mod . 2>/dev/null || true

$(OBJ)/read_input.o: read_input.f90 $(OBJ)/rjmcmc_com.o
	$(FC) -c -o $(OBJ)/read_input.o read_input.f90

$(OBJ)/alloc_obj.o: alloc_obj.f90 $(OBJ)/rjmcmc_com.o
	$(FC) -c -o $(OBJ)/alloc_obj.o alloc_obj.f90

$(OBJ)/loglhood.o: loglhood.f90 $(OBJ)/quicksort.o $(OBJ)/rjmcmc_com.o
	$(FC) -c -o $(OBJ)/loglhood.o loglhood.f90

$(OBJ)/UPDATE_MAPfile.o: UPDATE_MAPfile.f90 $(OBJ)/rjmcmc_com.o
	$(FC) -c -o $(OBJ)/UPDATE_MAPfile.o UPDATE_MAPfile.f90

$(OBJ)/xcorr.o: xcorr.f90 $(OBJ)/data_type.o
	$(FCSERIAL) -c -o $(OBJ)/xcorr.o xcorr.f90

$(OBJ)/COVmatEst.o: COVmatEst.f90 $(OBJ)/data_type.o
	$(FCSERIAL) -c -o $(OBJ)/COVmatEst.o COVmatEst.f90

$(OBJ)/UpdateCOV.o: UpdateCOV.f90 $(OBJ)/rjmcmc_com.o
	$(FC) -c -o $(OBJ)/UpdateCOV.o UpdateCOV.f90

$(OBJ)/nrtype.o: nrtype.f90
	$(FCSERIAL) -c -o $(OBJ)/nrtype.o nrtype.f90

$(OBJ)/nr.o: nr.f90 $(OBJ)/nrtype.o
	$(FCSERIAL) -c -o $(OBJ)/nr.o nr.f90

$(OBJ)/nrutil.o: nrutil.f90 $(OBJ)/nrtype.o
	$(FCSERIAL) -c -o $(OBJ)/nrutil.o nrutil.f90

$(OBJ)/quicksort.o: quicksort.f90
	$(FCSERIAL) -c -o $(OBJ)/quicksort.o quicksort.f90
