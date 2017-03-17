include Makefile.opts

SRC_DIR   = ../tools
SRC_DIR2  = tools
BUILD_DIR = build
DATA_DIR  = data
EXEC      = laplace
SOURCES   = laplace.c scalapackinfo.c LaplaceMatinit.c stringTools.c pblasIOtools.c nameFile.c
OBJECTS   = $(SOURCES:%.c=$(BUILD_DIR)/%.o)
LIB       = $(SCALAPACK)
INC       = $(INCBLAS) $(INCTOOLS)

COMPILE   = $(CC) $(CFLAGS) $(INC) -c
LINK      = $(F77) $(LDFLAGS)

all: $(OBJECTS) laplace

%.o: %.c
	$(COMPILE) $*.c -o $@

$(BUILD_DIR)/%.o: %.c
	$(COMPILE) $< -o $@
	
$(BUILD_DIR)/%.o: $(SRC_DIR2)/%.c
	$(COMPILE) $< -o $@

laplace: $(OBJECTS) 
	$(LINK) $(OBJECTS) $(LIB) -lm -o $@

clean:
	rm -f $(EXEC) $(BUILD_DIR)/*.o $(DATA_DIR)/*.out
