#############################################################################
## Makefile                                                                ##
##                                                                         ##
## Makefile for compiling using g++                                        ##
##                                                                         ##
##    11/26/04 WQ -- Makefile created                                      ##
##    12/09/19 WQ -- Makefile updated                                      ##
#############################################################################
PRJT   = jacobiBC_OMP_Test
NAME   = testJacobiBC_OMP
SRCT   = jacobiBC_OMP.c jacobiBC_OMP_Test.c
INCL   = jacobi.h const.h
OBJT   = $(SRCT:%.c=%.o)

CC     = g++
CFLAGS = -O3 -fopenmp
#LFLAGS = -lm

# the default target
$(PRJT) : $(OBJT)
	$(CC) -g $(CFLAGS) -o $(NAME) $(OBJT) -lm 

%.o : %.c $(INCL)
	$(CC) -c -g $(CFLAGS) -o $@ $<

# removes files created during build
clean:
	rm -f $(OBJT) $(PRJT)
