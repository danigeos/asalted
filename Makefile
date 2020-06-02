CC = gcc
BIN = ~/software/bin
LIBS	= -lm -lc 
C_OPTS	= -g

asalted:
%:	
	$(CC) $*.c -o $(BIN)/$@ $(C_OPTS) $(LIBS) 
	chmod  a+x  $(BIN)/$@
