##comments
ICC= g++
DEBUG=
CFLAGS= -O3 -Wall -m64 -g -Wno-deprecated

PNAME = stuff
LIBS = 
#-lpthread 
#

default : fort

fort: precompile $(PNAME)

%.o: %.cpp
	$(ICC) $(CFLAGS) -c $*.cpp

precompile : $(PRE)

GregsMathHdr.o : GregsMathHdr.cpp GregsMathHdr.h
	$(ICC) $(CFLAGS) -c GregsMathHdr.cpp
	
mcmchdr.o : mcmchdr.cpp mcmchdr.h random.h GregsMathHdr.h 
	$(ICC) $(CFLAGS) -c mcmchdr.cpp
	
all_start_spherical_jeans.o : all_start_spherical_jeans.cpp los_sph_jeans.h sph_jeans_new.h fun_sph_potential.h fun_beta.h fun_sph_stellar.h
	$(ICC) $(CFLAGS) -c all_start_spherical_jeans.cpp
	
ALL_SPH_OBJS_BOOST = mcmchdr.o  GregsMathHdr.o AndrewsMathHdr.o  all_start_spherical_jeans.o 
sph_boost:	$(ALL_SPH_OBJS_BOOST) 
	$(ICC) -I /usr/local/boost_1_55_0/  -o ~/Research/jvalue_project/node_stuff//sph_boost $(ALL_SPH_OBJS_BOOST) $(LIBS) $(CFLAGS)  \/usr/local/lib/libboost_program_options.a
	
NODE_TEST_OBJ = AndrewsMathHdr.o GregsMathHdr.o node_test.o
node_test:	$(NODE_TEST_OBJ) 
	$(ICC) -I /usr/local/boost_1_55_0/  -o ~/Research/jvalue_project/node_stuff/test_exe $(NODE_TEST_OBJ) $(LIBS) $(CFLAGS)  \/usr/local/lib/libboost_program_options.a


clean:
	rm -f *.o *.mod *.d *.pc *.obj *~ $(PNAME)

