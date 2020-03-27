CC = gcc
LD = gcc
CFLAGS = -Ofast -march=native
INCLUDES=
LDFLAGS =
LDLIBS = -lm
OBJS = barnes_hut.o
SRC = barnes_hut.c
TASK = barnes_hut
#all: $(TASK)


$(TASK): $(OBJS)$
	$(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(TASK)

run: $(TASK)$
	./$(TASK) 10 "input_data/ellipse_N_00010.gal" 200 0.00001 0.255 0

sun: $(TASK)$
	./$(TASK) 4 "input_data/sun_and_planets_N_4.gal" 10000 0.00001 0

test: $(TASK)$
	./$(TASK) 2000 "input_data/ellipse_N_02000.gal" 200 0.00001 0.255 0
	cd compare_files && ./a.out 2000 "../result.gal" "../ref_output_data/ellipse_N_02000_after200steps.gal" && cd ..


${OBJS}: $(SRC)
	$(CC) $(INCLUDES) $(CFLAGS) -c $(SRC) #-lm

clean:
	$(RM) $(TASK) $(OBJS)
