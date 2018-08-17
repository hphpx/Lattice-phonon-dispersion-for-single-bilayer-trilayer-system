TARGET=test
INCLUDES=-I /home/hong/eigenC++/eigen-eigen-b3f3d4950030/

$(TARGET): trilayertri.o DynTrilayer.o trimain.o
	g++ -o $@ trilayertri.o DynTrilayer.o trimain.o
trimain.o:   trimain.cpp
	g++ -c trimain.cpp   $(INCLUDES)
trilayertri.o: trilayertri.cpp
	g++ -c trilayertri.cpp $(INCLUDES)
DynTrilayer.o: DynTrilayer.cpp
	g++ -c DynTrilayer.cpp -lm
clean:
	rm *.o  $(TARGET)
