# Project: G-SOMAnts
# Makefile created by Dev-C++ 5.11

CPP      = g++
CC       = gcc
#WINDRES  = windres.exe
RES      =
#G-SOMAnts_private.res
OBJ      = main.o aloca.o functions.o $(RES)
LINKOBJ  = main.o aloca.o functions.o $(RES)
LIBS     =
#-L"C:/Program Files (x86)/Dev-Cpp/MinGW64/lib32" -L"C:/Program Files (x86)/Dev-Cpp/MinGW64/x86_64-w64-mingw32/lib32" -static-libgcc -lgdi32  -lcomdlg32  -luuid  -loleaut32  -lole32 -m32
INCS     =
#-I"C:/Program Files (x86)/Dev-Cpp/MinGW64/include" -I"C:/Program Files (x86)/Dev-Cpp/MinGW64/x86_64-w64-mingw32/include" -I"C:/Program Files (x86)/Dev-Cpp/MinGW64/lib/gcc/x86_64-w64-mingw32/4.9.2/include"
CXXINCS  =
#-I"C:/Program Files (x86)/Dev-Cpp/MinGW64/include" -I"C:/Program Files (x86)/Dev-Cpp/MinGW64/x86_64-w64-mingw32/include" -I"C:/Program Files (x86)/Dev-Cpp/MinGW64/lib/gcc/x86_64-w64-mingw32/4.9.2/include" -I"C:/Program Files (x86)/Dev-Cpp/MinGW64/lib/gcc/x86_64-w64-mingw32/4.9.2/include/c++"
BIN      = SSPSO
CXXFLAGS = $(CXXINCS)
# -m32
CFLAGS   = $(INCS)
# -m32
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

main.o: main.c
	$(CPP) -c main.c -o main.o $(CXXFLAGS)

aloca.o: aloca.c
	$(CPP) -c aloca.c -o aloca.o $(CXXFLAGS)

functions.o: functions.c
	$(CPP) -c functions.c -o functions.o $(CXXFLAGS)
