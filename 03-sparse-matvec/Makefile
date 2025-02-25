CC1 = gcc
CC2 = mpicc
CFLAGS = -Wall -O2
TARGET1 = mat-vec
TARGET2 = mat-vec-p
OBJ1 = mat-vec.o
OBJ2 = mat-vec-p.o

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJ1)
	$(CC1) $(CFLAGS) -o $(TARGET1) $(OBJ1)

$(TARGET2): $(OBJ2)
	$(CC2) $(CFLAGS) -o $(TARGET2) $(OBJ2)

mat-vec.o: mat-vec.c files.h
	$(CC1) $(CFLAGS) -c mat-vec.c

mat-vec-p.o: mat-vec-p.c files.h
	$(CC2) $(CFLAGS) -c mat-vec-p.c

clean:
	rm -f $(TARGET1) $(TARGET2) $(OBJ1) $(OBJ2) *.txt

.PHONY: all clean