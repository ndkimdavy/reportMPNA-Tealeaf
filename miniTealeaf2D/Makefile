# ===== miniTealeaf2D =====

FC     = gfortran
FFLAGS = -O2 -Wall -std=f2008
TARGET = miniTealeaf2D

OBJ = config.o hydro.o radiation.o io.o main.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ)

config.o: config.f90
	$(FC) $(FFLAGS) -c config.f90

hydro.o: hydro.f90 config.o
	$(FC) $(FFLAGS) -c hydro.f90

radiation.o: radiation.f90 config.o
	$(FC) $(FFLAGS) -c radiation.f90

io.o: io.f90 config.o
	$(FC) $(FFLAGS) -c io.f90

main.o: main.f90 config.o hydro.o radiation.o io.o
	$(FC) $(FFLAGS) -c main.f90

clean:
	rm -f *.o *.mod $(TARGET)