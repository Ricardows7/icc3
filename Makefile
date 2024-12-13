# Programas
PROGS = ajustaPol ajustaPolMelhorado gera_entrada

# Compilador e flags
CC     = gcc
CFLAGS = -Wall -g -O3 -mavx -march=native -DLIKWID_PERFMON -I${LIKWID_HOME}/include
LFLAGS = -lm -llikwid -L${LIKWID_LIB}

# Objetos
OBJ_AJUSTAPOL = ajustaPol.o utils.o
OBJ_AJUSTAPOLMELHORADO = ajustaPolMelhorado.o utils.o
OBJ_GERA = gera_entrada.o utils.o

all: $(PROGS)

# Regras para ajustaPol
ajustaPol: $(OBJ_AJUSTAPOL)
	$(CC) -o $@ $^ $(LFLAGS)

ajustaPol.o: ajustaPol.c utils.h
	$(CC) -c $< $(CFLAGS)

# Regras para ajustaPolMelhorado
ajustaPolMelhorado: $(OBJ_AJUSTAPOLMELHORADO)
	$(CC) -o $@ $^ $(LFLAGS)

ajustaPolMelhorado.o: ajustaPolMelhorado.c utils.h
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
