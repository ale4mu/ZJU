all: src/main.cpp src/mandelbrot.cpp src/window_show.cpp src/test.cpp
	cd src && g++ -fopenmp -o main main.cpp mandelbrot.cpp window_show.cpp `pkg-config --cflags --libs opencv4`
	cd src && g++ -o test test.cpp mandelbrot.cpp window_show.cpp `pkg-config --cflags --libs opencv4`

doc: doc/report.tex doc/mandelbrot.bib 
	doxygen Doxyfile
	cd doc && xelatex report.tex
	cd doc && bibtex report
	cd doc && xelatex report.tex
	cd doc && xelatex report.tex



clean:
	rm -f doc/*.aux doc/*.log doc/*.out doc/*.toc doc/*.bbl doc/*.blg doc/*.lof doc/*.lot doc/*.pdf 
	rm -rf doc/html
	rm -f src/main src/test src/*.png src/*.jpg src/*.bmp src/*.PNG src/*.JPG src/*.BMP src/*.mp4

.PHONY: all clean doc