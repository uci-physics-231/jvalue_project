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

NODE_TEST = node_test.o
sph_boost:	$(ALL_SPH_OBJS_BOOST) 
	$(ICC) -I /usr/local/boost_1_55_0/  -o ~/Research/jvalue_project/node_stuff/jeans $(NODE_TEST) $(LIBS) $(CFLAGS)  \/usr/local/lib/libboost_program_options.a


clean:
	rm -f *.o *.mod *.d *.pc *.obj *~ $(PNAME)

