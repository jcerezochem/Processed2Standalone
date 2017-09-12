SHELL=/bin/bash
INSTALL_DIR=/home/cerezo/mis_bin/
PROGRAM=procesed2stadalonetop.py

all: update_git_flag copy_to_installdir

update_git_flag:
	@git_hash=$$(git describe --long --dirty --always); \
	 git_date=$$(git show -s --format=%ci); \
	 git_hash_old=$$(grep "COMMIT = " $(PROGRAM)); \
	 git_date_old=$$(grep "DATE = " $(PROGRAM)); \
	 echo "COMMIT: $$git_hash vs $$git_hash_old)"

copy_to_installdir:
	cp $(PROGRAM) $(INSTALL_DIR) 

