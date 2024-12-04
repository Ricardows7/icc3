# Programas
PROGS = ajustePol gera_entrada

# Compilador e flags
CC     = gcc
CFLAGS = -Wall -g -O3 -mavx -march=native -DLIKWID_PERFMON -I${LIKWID_HOME}/include
LFLAGS = -lm -llikwid -L${LIKWID_LIB}

# Objetos
OBJ_AJUSTE = ajustePol.o utils.o
OBJ_GERA = gera_entrada.o utils.o

all: $(PROGS)

# Regras para ajustePol
ajustePol: $(OBJ_AJUSTE)
	$(CC) -o $@ $^ $(LFLAGS)

ajustePol.o: ajustePol.c utils.h
	$(CC) -c $< $(CFLAGS)

# Regras para gera_entrada
gera_entrada: $(OBJ_GERA)
	$(CC) -o $@ $^ $(LFLAGS)

gera_entrada.o: gera_entrada.c
	$(CC) -c $< $(CFLAGS)

# Regra para utils.c
utils.o: utils.c utils.h
	$(CC) -c $< $(CFLAGS)

# Limpeza
clean:
	@echo "Limpando arquivos objetos e executáveis..."
	@rm -f *.o $(PROGS)

purge: clean
	@echo "Limpando arquivos adicionais..."
	@rm -f *~ *.bak core

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tar)..."
	@ln -s . $(DISTDIR)
	@tar -cvf $(DISTDIR).tar $(addprefix ./$(DISTDIR)/, $(DISTFILES))
	@rm -f $(DISTDIR)

