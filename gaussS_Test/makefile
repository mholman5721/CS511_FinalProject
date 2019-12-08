#############################################################################
## Makefile                                                                ##
##                                                                         ##
## Makefile for compiling using g++                                        ##
##                                                                         ##
##    11/26/04 WQ -- Makefile created                                      ##
##    12/04/19 WQ -- Makefile updated                                      ##
#############################################################################
PRJT   = gaussS_Test
NAME   = testGaussS
SRCT   = gaussS.c gaussS_Test.c
INCL   = gaussS.h const.h
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
