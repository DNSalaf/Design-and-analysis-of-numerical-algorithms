CPP = $(wildcard src/*.c)
EXE = output

all:
	gcc $(CPP) -o $(EXE)

run:
	./$(EXE)

clean:
	rm -f $(EXE)