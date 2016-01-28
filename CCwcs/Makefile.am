## Process this file with automake to produce Makefile.in
#
# Author: Josh Cardenzana, jvcard@iastate.edu
# Date: Oct, 2015
#

AM_CXXFLAGS = -std=c++11 -I./include $(WCS_CXXFLAGS)
AM_CFLAGS = 

lib_LTLIBRARIES = libccwcs.la
libccwcs_la_SOURCES = WcsHandler.cxx
libccwcs_la_HEADERS = WcsHandler.h
libccwcs_la_CPPFLAGS =  $(CCWCS_CPPFLAGS)
libccwcs_la_LIBADD=$(LIBS)
libccwcs_ladir = $(libdir)
libccwcs_la_DATA = WcsHandler.cxx

WcsHandler.cxx: $(libwcs_la_HEADERS)

#bin_PROGRAMS = testwcs
#testwcs_SOURCES = main.cpp
#testwcs_LDADD = libccwcs.la
#testwcs_CPPFLAGS = $(CCWCS_CPPFLAGS)

CLEANFILES = 