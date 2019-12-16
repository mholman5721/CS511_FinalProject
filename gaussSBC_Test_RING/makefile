#############################################################################
## Makefile                                                                ##
##                                                                         ##
## Makefile for compiling using g++                                        ##
##                                                                         ##
##    11/26/04 WQ -- Makefile created                                      ##
##    12/09/19 WQ -- Makefile updated                                      ##
#############################################################################
PRJT   = testGaussSBC_Ring
NAME   = testGaussSBC_Ring
SRCT   = gaussSBC.c gaussSBC_Test.c
INCL   = gaussSBC.h const.h
OBJT   = $(SRCT:%.c=%.o)

CC     = g++
CFLAGS = -O3
#LFLAGS = -lm

# the default target
$(PRJT) : $(OBJT)
	$(CC) $(CFLAGS) -o $(NAME) $(OBJT)

%.o : %.c $(INCL)
	$(CC) -c $(CFLAGS) -o $@ $<

# removes files created during build
clean:
	rm -f $(OBJT) $(PRJT)
