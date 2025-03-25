# Compiler and flags
CC       = mpicc
CFLAGS   = -O2 -Wall -I/usr/local/include
LDFLAGS  = -L/usr/local/lib -lHYPRE -lm

# Source files and objects
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SOURCES  = $(wildcard $(SRCDIR)/*.c)
OBJECTS  = $(patsubst $(SRCDIR)/%.c, $(OBJDIR)/%.o, $(SOURCES))
TARGET   = $(BINDIR)/main

# Create necessary directories
$(shell mkdir -p $(OBJDIR) $(BINDIR))

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Object file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -rf $(OBJDIR)/*.o $(TARGET)

.PHONY: all clean
