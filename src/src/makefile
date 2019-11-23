#############################################################################
## Makefile                                                                ##
##                                                                         ##
## Makefile for compiling on Linux/Unix systems using gcc/cc               ##
##                                                                         ##
##    11/26/04 WQ -- Makefile created                                      ##
#############################################################################
PRJT   = jacobiBC
SRCT   = jacobi.c jacobiBC.c
INCL   = jacobi.h const.h
OBJT   = $(SRCT:%.c=%.o)

CC     = cc
CFLAGS = -O2
LFLAGS = -lm

# the default target
$(PRJT) : $(OBJT)
	$(CC) $(LFLAGS) -o $(PRJT) $(OBJT)

%.o : %.c $(INCL)
	$(CC) -c $(CFLAGS) -o $@ $<

# removes files created during build
clean:
	rm -f $(OBJT) $(PRJT)
