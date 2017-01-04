TPROGS = getTargetsDef getAccssnTaxID getfilesToTaxNodes #getGInTaxID
PROGS = cuCLARK cuCLARK-l $(TPROGS)

.PHONY: all clean target_definition

# install all programs in folder ./exe/
all:
	$(MAKE) -C src
	@mkdir -p exe
	@cp $(addprefix src/,$(PROGS)) exe/

clean:
	rm -rf exe
	$(MAKE) -C src clean

target_definition:
	$(MAKE) -C src target_definition
	@mkdir -p exe
	@cp  $(addprefix src/,$(TPROGS)) exe/
