PROJ_NAME = project1
CC = gcc #arm-linux-gnueabi-gcc
VECTFLAGS = -ftree-vectorize -ffast-math -fsingle-precision-constant -mvectorize-with-neon-quad -ftree-vectorizer-verbose=6
CFLAGS = -Wall -g -O3 -march=armv7-a -mcpu=cortex-a8  -mfloat-abi=softfp -mfpu=neon $(VECTFLAGS) #-funroll-loops 
LIBS = -lm -lrt
OBJFILES := $(patsubst %.c,%.o,$(wildcard *.c))
DFLAG := -DINTRINSICS
ifeq ($(ASM), 1)
	OBJFILES += $(patsubst %.s, %.o, $(wildcard *.s))
	DFLAG = -DASM
endif
$(PROJ_NAME): $(OBJFILES) 
#	echo $(OBJFILES)
	$(CC) -o $(PROJ_NAME) $(OBJFILES) $(LIBS)
%.o: %.c
	$(CC) $(CFLAGS) $(DFLAG) -c -o $@ $<
%.o: %.s
	$(CC) $(CFLAGS) -c -o $@ $<
%.lst: %.c
	$(CC) -g $(CFLAGS) -Wa,-adhln $(LIBS) $< > $@
clean:
	rm -f *.o *.lst
