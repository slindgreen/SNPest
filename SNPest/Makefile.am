## Source directory

bin_PROGRAMS = EvoFoldV2 grammarTrain dfgEval dfgTrain multinomial dfgEval_SNPest

LDADD = $(top_srcdir)/phy/libphy.la -lboost_program_options -llapack -lntl -lopt -lnewmat -lm

EvoFoldV2_SOURCES    = EvoFold.cpp
grammarTrain_SOURCES = grammarTrain.cpp
dfgEval_SOURCES      = dfgEval.cpp
dfgTrain_SOURCES     = dfgTrain.cpp
multinomial_SOURCES  = multinomial.cpp
dfgEval_SNPest_SOURCES = dfgEval_SNPest.cpp

#  compiler options
AM_CPPFLAGS = -I$(top_srcdir)