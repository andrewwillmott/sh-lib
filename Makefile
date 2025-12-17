#
# Makefile for SH library
#

ifeq ($(origin CXX), default) # gmake defaults this to g++ instead of c++
	CXX = c++
endif
CXXFLAGS ?= --std=c++11
OPTS=-O3 -Wall
DBG_OPTS=-DVL_DEBUG -g

PREFIX ?= /usr/local
LIB_DIR = $(PREFIX)/lib
INCLUDE_DIR = $(PREFIX)/include

SHL_HEADERS = SHLib.hpp ZHLib.hpp
SHL_SOURCES = SHLib.cpp ZHLib.cpp
SHL_DEPS    = $(SHL_HEADERS) VL234f.hpp Makefile
SHL_OBJS    = $(SHL_SOURCES:.cpp=.o)
SHL_DOBJS   = $(SHL_SOURCES:.cpp=D.o)

all: libsh.a libshd.a

test: libshd.a ${SHL_HEADERS} SHLibTest.cpp
	$(CXX) $(CXXFLAGS) $(DBG_OPTS) -o $@ -L. -lshd SHLibTest.cpp
	./$@

libsh.a: $(SHL_OBJS)
	$(AR) rcs $@ $^

libshd.a: $(SHL_DOBJS)
	$(AR) rcs $@ $^

%.o: %.cpp $(SHL_HEADERS)
	$(CXX) $(CXXFLAGS) $(OPTS) -c $< -o $@

%D.o: %.cpp $(SHL_HEADERS)
	$(CXX) $(CXXFLAGS) $(DBG_OPTS) -c $< -o $@

install: libsh.a libshd.a
	mkdir -p $(LIB_DIR)
	cp libsh.a libshd.a $(LIB_DIR)
	mkdir -p $(INCLUDE_DIR)
	cp $(SHL_HEADERS) $(INCLUDE_DIR)

clean:
	rm -rf *.o *.a *.dSYM test
