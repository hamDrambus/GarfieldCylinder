INCDIR = $(GARFIELD_HOME)/include
HEEDDIR = $(GARFIELD_HOME)/share/Heed
LIBDIR = $(GARFIELD_HOME)/lib
TARNAME = CylinderCharge

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
#CFLAGS += -g

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` `root-config --ldflags` -lGeom \
	-lgfortran -lm

#LDFLAGS += -g

target: $(TARNAME).o MaterialCustomDrift.o
	$(CXX) $(TARNAME).o MaterialCustomDrift.o -o $(TARNAME) $(LDFLAGS)
	rm -f $(TARNAME).o MaterialCustomDrift.o

$(TARNAME).o: $(TARNAME).cpp
	$(CXX) $(CFLAGS) $(TARNAME).cpp
	
MaterialCustomDrift.o: MaterialCustomDrift.cpp MaterialCustomDrift.hh
	$(CXX) $(CFLAGS) MaterialCustomDrift.cpp

clean:
	rm $(TARNAME).o MaterialCustomDrift.o $(TARNAME)
