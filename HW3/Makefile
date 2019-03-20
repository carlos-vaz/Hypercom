CC = gcc
LFLAGS = -lpthread
OBJ = main

build = \
	@echo "[BUILDING $@...]"; \
	$1; 

% : %.c
	$(call build, $(CC) $< -o $@ $(LFLAGS))

all: $(OBJ)

clean:
	rm $(OBJ) 

