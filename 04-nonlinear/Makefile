CC = mpicc
CFLAGS = -Wall -O0 -g3 
LDFLAGS = -lhypre -lm
SRC = flame.c
OBJ = flame.o
EXEC = flame

$(EXEC): $(OBJ)
    $(CC) -o $(EXEC) $(OBJ) $(LDFLAGS)

$(OBJ): $(SRC)
    $(CC) -c $(SRC) $(CFLAGS)

clean:
    rm -f $(OBJ) $(EXEC)
