#
# This makefile is divided into four sections:
#    (1) objects
#    (2) sources
#    (3) headers
#    (4) subdirectories
#
# Instructions for how to define these is given in each section below
#
include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

##################### (1) objects
#
# (1.1) If there are objects in this directory:
#    (1) define the OBJ macro as the list of object files
#    (2) uncomment the OBJ_LINK and DEP macros here
#    (3) comment the OBJ, OBJ_LINK, and DEP macros in the
#        section (1.1).
#
OBJ = \
	FSDEMatSupport2DT.o \
	FSDEMat2DT.o \
	Mech_PK2_2D.o \
	Mech_TanMod_2D.o \
	ME_MixedModulus_2D.o \
	ME_MixedModulus_2Dspatial.o \
	ME_PK2_2D.o \
	ME_TanMod_2D.o \
	Elec_PK2_2D.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)

#
# (1.2) If there are not object in this directory:
#    (1) comment the OBJ, OBJ_LINK, and DEP definitions in section (1.1)
#    (2) uncomment the OBJ, OBJ_LINK, and DEP definitions below
#
#OBJ = dummy
#OBJ_LINK =
#DEP = $(MACRO_DIR)/empty

##################### (2) sources
#
# (2) If there are sources in this directory list them in the
#     SRC macro; otherwise, comment the SRC macro.
#
SRC = \
	FSDEMatSupport2DT.cpp \
	FSDEMat2DT.cpp \
	Mech_PK2_2D.cpp \
	Mech_TanMod_2D.cpp \
	ME_MixedModulus_2D.cpp \
	ME_MixedModulus_2Dspatial.cpp \
	ME_PK2_2D.cpp \
	ME_TanMod_2D.cpp \
	Elec_PK2_2D.cpp

##################### (3) headers
#
# (3) If there are header files in this directory:
#    (1) list the header files in the HDR macro
#    (2) uncomment the HDR_LINK
# 
# If there are no header files in this directory comment both the
# HDR and HDR_LINK macro   
#
HDR = \
	FSDEMatSupport2DT.h \
	FSDEMatSupport2DT.i.h \
	FSDEMat2DT.h \
	FSDEMat2DT.i.h \
	FSDE_common_defines2D.h \
	FSDE_inc2D.h 

HDR_LINK = $(HDR:.h=.h_link)

##################### (4) sub-directories
#
# (4.1) If there are sub-directories contained in this directory:
#    (1) list the subdirectories in the SUB_DIR macro
#    (2) uncomment the subdir_driver macro
#    (3) comment the SUB_DIR and subdir_driver in section (4.2)
#
#SUB_DIR = 
#subdir_driver: subdir_loop

#
# (4.2) If there are no sub-directories in this directory:
#    (1) comment the SUB_DIR and subdir_driver in section (4.1)
#    (2) uncomment the subdir_driver macro below
#
SUB_DIR =
subdir_driver: 

#####################

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = $(MACRO_DIR)/empty
include $(DEPEND)
