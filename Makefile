#---------------------------------------- Special behavior(all closed as default)
#OPT += -DTURN_OFF_MFT                    # turn off ddm mf transfer
#OPT += -DTURN_OFF_CMT                    # turn off ddm cm transfer
OPT += -DVERBOSE                          # More Terminal output
OPT += -DMAC                              # For Mac system

#---------------------------------------- Modes
OPT += -DEXPLORE_PARAM_SPACE
#OPT += -DOUTPUT_TK
#--------------------------------------- Select target computer
#SYSTYPE="dlcheng"
SYSTYPE="Mac"
#--------------------------------------- Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl-1.16/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl-1.16/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib  -Wl
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   =  ddm_tk

OBJS   =  allvars.o bias_mf_norm.o cdm_bias.o cdm_mf.o cl.o cm_relation.o ddm_bias.o ddm_mf.o ddm_transfer.o \
          decay_param.o filters.o growth.o init_all.o main.o mass_mapping.o nfw_profile.o power_cal.o power_ui.o\
          set_params.o smith2.o smooth_field.o spline_2d.o tk.o tk_bbks.o tools.o variance.o

INCL   = allvars.h proto.h define.h smith2.h spline_2d.h Makefile

CFLAGS = $(OPTIONS) $(GSL_INCL) 

LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC) *.gch
