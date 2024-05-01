SHELL=/bin/bash
INSTALL_DIR=/home/cerezo/mis_bin/
PROGRAM=procesed2stadalonetop.py
VERSION_FILE=version_tag.txt

all: update_git_flag copy_to_installdir

update_git_flag:
	@echo "Updating version flag: $(VERSION_FILE)"
	@if [ ! -e $(VERSION_FILE) ]; then                   \
	     echo "class version_tag:" > $(VERSION_FILE);    \
	     echo "    COMMIT = \"XXX\""   >> $(VERSION_FILE);    \
	     echo "    DATE = \"000\""    >> $(VERSION_FILE);    \
	 fi
	@git_hash=$$(git describe --long --dirty --always); \
	 git_hash_old=$$(grep "COMMIT = " $(VERSION_FILE)); \
	 git_hash_old=$${git_hash_old##*COMMIT = };\
	 git_hash_old=$${git_hash_old//\"/};\
	 git_date=$$(git show -s --format=%ci); \
	 git_date_old=$$(grep "DATE = " $(VERSION_FILE)); \
	 git_date_old=$${git_date_old##*DATE = };\
	 git_date_old=$${git_date_old//\"/};\
	 if [ "$$git_hash" != "$$git_hash_old" ]; then \
	     echo "Git commit has changed:"; \
	     echo "COMMIT: $$git_hash - OLD: $$git_hash_old"; \
	     echo "DATE:   $$git_date - OLD: $$git_date_old"; \
	     sed -i "s/$$git_hash_old/$$git_hash/" $(VERSION_FILE); \
	     sed -i "s/$$git_date_old/$$git_date/" $(VERSION_FILE); \
	 fi

copy_to_installdir:
	@echo "Installing to: $(INSTALL_DIR)"
	@cp $(PROGRAM) $(INSTALL_DIR) 
	@sed -i "9,11d" $(INSTALL_DIR)/$(PROGRAM)
	@sed -i "8r$(VERSION_FILE)" $(INSTALL_DIR)/$(PROGRAM)

