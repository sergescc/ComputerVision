PART1 = HW3a

PART2 = HW3b

CFLAGS = -g $(INCLUDE) -Wall -lm

INCLUDE = -I inc

LIBS = -pthread

CC = gcc

SOURCES1 = Homework3a.c CursorCntl.c ImageProcessing.c

SOURCES2 = Homework3b.c CursorCntl.c ImageProcessing.c

TARGETDIR = bin

MODULES = $(TARGETDIR)/modules

OBJS1 = $(addprefix $(MODULES)/,$(SOURCES1:.c=.o))

OBJS2 = $(addprefix $(MODULES)/,$(SOURCES2:.c=.o))

vpath %.c ./src/
vpath %.h ./inc/

all: $(TARGETDIR)/$(PART1) $(TARGETDIR)/$(PART2)

$(TARGETDIR)/$(PART1): $(OBJS1)
	if [ ! -d $(TARGETDIR) ]; then \
		mkdir -p $(TARGETDIR); \
	fi
	$(CC) $(CFLAGS) $(LIBS) $(OBJS1) -o $@

$(TARGETDIR)/$(PART2): $(OBJS2)
	if [ ! -d $(TARGETDIR) ]; then \
		mkdir -p $(TARGETDIR); \
	fi
	$(CC) $(CFLAGS) $(LIBS) $(OBJS2) -o $@

$(MODULES)/%.o: %.c
	if [ ! -d $(MODULES) ]; then \
		mkdir -p $(MODULES); \
	fi
	$(CC) $(CFLAGS) $(LIBS) -c -o $@ $<


clean:
	rm $(OBJS) $(TARGETDIR)/$(PART1) $(TARGETDIR)/$(PART2)
