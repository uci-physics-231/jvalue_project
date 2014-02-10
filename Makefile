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

NODE_TEST_OBJ = AndrewsMathHdr.o GregsMathHdr.o node_test.o
node_test:	$(NODE_TEST_OBJ) 
	$(ICC) -I /usr/local/boost_1_55_0/  -o ~/Research/jvalue_project/node_stuff/test_exe $(NODE_TEST_OBJ) $(LIBS) $(CFLAGS)  \/usr/local/lib/libboost_program_options.a


clean:
	rm -f *.o *.mod *.d *.pc *.obj *~ $(PNAME)

