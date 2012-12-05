include $(GWOLLUM_LOCATION)/src/include.mk
OINSTALL_LIB   = ../lib
OINSTALL_BIN   = ../bin

Omicron_LIB       = libOmicron.so
Omicron_OBJS      = Omicron.o OmicronDict.o
Omicron_DICT_HDRS = Omicron.h

omicron_EXE       = omicron.exe
omicron_OBJS      = omicron.o

CXXFLAGS += -I$(UTILS_SRC)
CXXFLAGS += -I$(SEGMENTS_SRC)
CXXFLAGS += -I$(TRIGGERS_SRC)
CXXFLAGS += -I$(TIMESLIDES_SRC)
CXXFLAGS += -I$(FRAMELIB_INCLUDE)
CXXFLAGS += -I$(FFTW3_INCLUDE)
CXXFLAGS += -I$(OMICRON_SRC)

all: $(Omicron_LIB) $(omicron_EXE)

%.o : %.cc
	$(CXX) -c $< $(CXXFLAGS)

$(Omicron_LIB) : $(Omicron_OBJS)
	$(CXX) -shared $^ -o $@ $(CXXFLAGS) $(ROOT_LIBS) $(FFTW3_LIB) -L$(FRAMELIB_LIB) -lFrame -L$(INSTALL_LIB) -lCUtils -lRootUtils -lTimeSlides -lSegments -lTriggers -lOmicronToolBox
	cp $@ $(OINSTALL_LIB)

OmicronDict.cc OmicronDict.h : $(Omicron_DICT_HDRS)
	$(ROOTSYS)/bin/rootcint -f OmicronDict.cc -c $(ROOT_INCLUDE) -I$(FFTW3_INCLUDE) -I$(FRAMELIB_INCLUDE) -I$(UTILS_SRC) -I$(SEGMENTS_SRC) -I$(TRIGGERS_SRC) -I$(TIMESLIDES_SRC) -I$(OMICRON_SRC) -p $(Omicron_DICT_HDRS) LinkDef.h

$(omicron_EXE) : $(omicron_OBJS)
	$(CXX) $^  -g -o $@ $(CXXFLAGS) $(ROOT_LIBS) $(FFTW3_LIB) -L$(FRAMELIB_LIB) -lFrame -L$(INSTALL_LIB) -lCUtils -lRootUtils -lTimeSlides -lSegments -lTriggers -lOmicronToolBox -L$(OINSTALL_LIB) -lOmicron
	cp $@ $(OINSTALL_BIN)

clean :
	- rm -f *Dict.* *.exe *.o *.so 

