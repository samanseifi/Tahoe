# $Id: makefile,v 1.5 2005/08/12 05:27:13 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	LU_serial_driver.o \
	LU_serial_driver_init.o \
	LU_serial_driver_factorize.o \
	LU_serial_driver_solve.o \
	LU_serial_driver_free.o	

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

# sources
SRC = \
	LU_serial_driver.c \
	LU_serial_driver_init.c \
	LU_serial_driver_factorize.c \
	LU_serial_driver_solve.c \
	LU_serial_driver_free.c	

# headers
HDR = LU_serial_driver.h LU_serial_driver_int.h
HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
