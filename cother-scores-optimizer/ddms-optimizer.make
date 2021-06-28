#############################################################################
# Makefile for building: makeddms2s++
# Scope:  global
# Version: 1.0.0
#############################################################################

DLEXT = so
TARGETNAME = ddms
TARGET = $(TARGETNAME).$(DLEXT)
LN = ln -fs
RM = rm

all: ddms-lib $(TARGET)

clean: clean-ddms-lib clean-target

ddms-lib:
	swig -c++ -perl ddms.i && \
	perl Makefile.PL && \
	$(MAKE)

$(TARGET): blib/arch/auto/ddms/$(TARGET)
	$(LN) $< $@

clean-ddms-lib:
	- $(MAKE) clean

clean-target:
	- $(RM) $(TARGET)

.PHONY: all clean ddms-lib clean-ddms-lib clean-target

