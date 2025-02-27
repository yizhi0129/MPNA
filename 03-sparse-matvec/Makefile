CC1 = gcc
CC2 = mpicc
CFLAGS = -Wall -O0 -g3

TARGET1 = mat-vec
TARGET2 = mat-vec-p
TARGET3 = eigenvalue
TARGET4 = eigenvalue-p
TARGET5 = gershgorin
TARGET6 = accuracy

OBJ1 = mat-vec.o
OBJ2 = mat-vec-p.o
OBJ3 = eigenvalue.o
OBJ4 = eigenvalue-p.o
OBJ5 = gershgorin.o
OBJ6 = accuracy.o

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6)

$(TARGET1): $(OBJ1)
	$(CC1) $(CFLAGS) -o $(TARGET1) $(OBJ1) -lm

$(TARGET2): $(OBJ2)
	$(CC2) $(CFLAGS) -o $(TARGET2) $(OBJ2) -lm

$(TARGET3): $(OBJ3)
	$(CC1) $(CFLAGS) -o $(TARGET3) $(OBJ3) -lm

$(TARGET4): $(OBJ4)
	$(CC2) $(CFLAGS) -o $(TARGET4) $(OBJ4) -lm

$(TARGET5): $(OBJ5)
	$(CC1) $(CFLAGS) -o $(TARGET5) $(OBJ5) -lm

$(TARGET6): $(OBJ6)
	$(CC1) $(CFLAGS) -o $(TARGET6) $(OBJ6) -lm

mat-vec.o: mat-vec.c files.h
	$(CC1) $(CFLAGS) -c mat-vec.c

mat-vec-p.o: mat-vec-p.c files.h
	$(CC2) $(CFLAGS) -c mat-vec-p.c

eigenvalue.o: eigenvalue.c files.h
	$(CC1) $(CFLAGS) -c eigenvalue.c

eigenvalue-p.o: eigenvalue-p.c files.h
	$(CC2) $(CFLAGS) -c eigenvalue-p.c

gershgorin.o: gershgorin.c files.h
	$(CC1) $(CFLAGS) -c gershgorin.c

accuracy.o: accuracy.c files.h
	$(CC1) $(CFLAGS) -c accuracy.c

clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6)

.PHONY: all clean