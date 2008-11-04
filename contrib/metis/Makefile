

default:
	(cd Lib ; make ) 
	(cd Programs ; make ) 
	(cd Test ; make ) 
	../../admin/mkinstalldirs ${prefix}/include/metis/
	../../admin/mkinstalldirs ${prefix}/lib
	for header in Lib/*.h; do install -m 644 $$header ${prefix}/include/metis/; done
	install -m 644 libmetis.a ${prefix}/lib

clean:
	(cd Lib ; make clean )
	(cd Programs ; make clean )
	(cd Test ; make clean )

realclean:
	(cd Lib ; make realclean )
	(cd Programs ; make realclean )
	(cd Test ; make realclean )
