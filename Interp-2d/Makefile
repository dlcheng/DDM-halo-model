SYSTYPE="dlcheng"
#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/Users/dalongcheng/Install/gsl-1.16/include
GSL_LIBS =  -L/Users/dalongcheng/Install/gsl-1.16/lib
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/users/student/s0902248/Lib/gsl-1.9/lib
GSL_INCL =  -I/users/student/s0902248/Lib/gsl-1.9/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   =  Spline_2d_test

OBJS   =  interp_2d_test.o Spline_2d.o 

INCL   =  Spline_2d.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) 


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch
