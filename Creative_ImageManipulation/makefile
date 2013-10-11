CC = g++
PROG = creative
SRCS = mainer.cpp
ifeq ($(shell uname),Darwin)
CFLAGS = -Wall -O2 -w -I . -I/opt/local/include/OpenEXR
LIBS = -L/opt/local/lib -l IlmImf -l Imath -l Half
else
CFLAGS = -Wall -g -w -I . -I/usr/include/OpenEXR
LIBS = -L/usr/local/lib -lIlmImf -lImath -lHalf
endif
all: $(PROG)
$(PROG): $(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)
clean:
	rm -f $(PROG)
