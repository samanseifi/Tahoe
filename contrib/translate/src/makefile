# $Id: makefile,v 1.10 2004/11/16 01:04:35 paklein Exp $
# object directory makefile - compilation instructions

include $(MACRO_DIR)/suffix
include $(MACRO_DIR)/$(ARCH).macros

# objects
OBJ = \
	TranslateIOManager.o \
	ExtractIOManager.o \
	ExtractNode.o \
	ExtractElement.o \
	ExtractQuad.o \
	PointPlots.o \
	MergeResults.o \
	JoinResults.o \
	Scroller.o \
	main.o

OBJ_LINK = $(OBJ)
DEP = $(OBJ:.o=.d)
#OBJ = dummy
#OBJ_LINK =
#DEP = /dev/null

# sources
SRC = \
	TranslateIOManager.cpp \
	ExtractIOManager.cpp \
	ExtractNode.cpp \
	ExtractElement.cpp \
	ExtractQuad.cpp \
	PointPlots.cpp \
	MergeResults.cpp \
	JoinResults.cpp \
	Scroller.cpp \
	main.cpp

# headers
HDR = \
	TranslateIOManager.h \
	ExtractIOManager.h \
	ExtractNode.h \
	ExtractElement.h \
	ExtractQuad.h \
	MergeResults.h \
	JoinResults.h \
	Scroller.h \
	PointPlots.h

HDR_LINK = $(HDR:.h=.h_link)

# subdirectories
SUB_DIR = 
#subdir_driver: subdir_loop
# SUB_DIR is empty
subdir_driver: 

# instructions for subdirectories
include $(MACRO_DIR)/subdir.targets

# dependencies
DEPEND = /dev/null
include $(DEPEND)
