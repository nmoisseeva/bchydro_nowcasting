FC = r.compile
CC = cc

FLAGS = -src

LIBS = ${ARMNLIB}/lib/Linux_pgi611/librmn.a

OBJS1 = \
main.o \
sf_gen.o \
get_inputs.o \
getninj_Z.o \
get_ME_GP.o \
sub_rotation.o \
llacar.o \
read_ME_FST_G.o \
cal_dist_sphere.o \
test_coords.o \
read_sfcstation.o \
source_pt.o \
calcoord_stn.o \
gc2rc.o \
rc2rp.o \
find_neighb.o \
findswpt.o \
MD_prep.o  \
init_para.o \
MD_sph.o \
findlc.o \
cal_cor.o

F90DECKS = \
main.f90 \
sf_gen.f90 \
get_inputs.f90 \
getninj_Z.f90 \
get_ME_GP.f90 \
sub_rotation.f90 \
llacar.f \
read_ME_FST_G.f90 \
cal_dist_sphere.f90 \
test_coords.f90 \
read_sfcstation.f90 \
source_pt.f90 \
calcoord_stn.f90 \
gc2rc.f90 \
rc2rp.f90  \
find_neighb.f90 \
findswpt.f90 \
MD_prep.f90 \
init_para.f90 \
MD_sph.f90 \
findlc.f90 \
cal_cor.f90

PROGRAM1 = main

default: $(PROGRAM1) 

all:  $(PROGRAM1) 

$(PROGRAM1): 
	$(FC) $(FLAGS) $(F90DECKS)  -librmn -o $(PROGRAM1)

clean:
	rm *.o; rm $(PROGRAM1) 
