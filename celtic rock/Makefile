CXXFLAGS = -Wall -g -std=c++11 -lm
CXXLIBS = -lopengl32 -lfreeglut
SOURCE = main.cpp
DFLAGS =
IFLAGS = 
LFLAGS =
ifeq ($(OS), Windows_NT)
	IFLAGS = -I"E:\Compilers\MSYS\include"
	LFLAGS = -L"E:\Compilers\MSYS\lib"
endif


OBJS = $(SOURCE:.cpp=.o)
DEPS = $(SOURCE:.cpp=.d)

release: $(OBJS) $(DEPS)
	$(CXX) $(OBJS) -o main.exe $(CXXLIBS)  $(LFLAGS)

%.o: %.cpp %.d
	$(CXX) $(DFLAGS) -c -MD $< -o $@ $(CXXFLAGS) $(IFLAGS)

%.d: %.cpp
	$(CXX) $(DFLAGS) -c -MD $<  $(CXXFLAGS) $(IFLAGS)

clean:
	$(RM) $(OBJS) $(DEPS)

-include $(DEPS)