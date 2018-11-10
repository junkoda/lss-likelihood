all: install

build: $(SRC)
	python3 setup.py config build_ext --inplace

install:
	python3 setup.py install


.PHONY: clean check

clean:
	rm -rf build lss_likelihood/*.so
