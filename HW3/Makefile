CC = gcc
CFLAGS = -std=c99
LFLAGS = -lpthread \
	 -lm
OBJ = main

build = \
	@echo "[BUILDING $@...]"; \
	$1; 

% : %.c
	$(call build, $(CC) $< -o $@ $(CFLAGS) $(LFLAGS))

all: $(OBJ)

clean:
	rm $(OBJ) 

