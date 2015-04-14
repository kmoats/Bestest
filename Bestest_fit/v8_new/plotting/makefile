NAME	= higgsanalysis.exe

OBJECTS = $(NAME:.exe=.o) plotting.o

CPPFLAGS = -I/usr/local/include -I$(shell root-config --incdir)

ROOTLIBS = $(shell root-config --libs)

HEADERS = style.h plotting.h

all:	$(NAME)

$(OBJECTS): $(HEADERS)

$(NAME): $(OBJECTS)
	g++ -o $(NAME) -O $(OBJECTS) $(ROOTLIBS) -L/usr/local/lib

.cpp.o:
	g++ -c -O $(CPPFLAGS) -I$(ROOTSYS)/include $*.cpp