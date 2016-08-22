

all:	clean
	cd sage_output/src ; make
	cd magma_output/src ; make
	cd vec_lzz_p_extra/src ; make
	cd lzz_pX_middle_product/src ; make
	cd lzz_pX_CRT/src ; make
	cd ZZ_CRT/src ; make
	cd mat_lzz_pX_extra/src ; make
	cd bivariate/src ; make
	cd cauchy_like_geometric_special/src ; make
	cd mosaic_hankel/src ; make
	cd ZZ_pX_CRT/src ; make

clean:
	rm -f lib/*
	rm -f include/*
	cd sage_output/src ; make clean
	cd magma_output/src ; make clean
	cd vec_lzz_p_extra/src ; make clean
	cd lzz_pX_middle_product/src ; make clean
	cd lzz_pX_CRT/src ; make clean
	cd ZZ_CRT/src ; make clean
	cd ZZ_pX_CRT/src ; make clean
	cd mat_lzz_pX_extra/src ; make clean
	cd bivariate/src ; make clean
	cd cauchy_like_geometric_special/src ; make clean
	cd mosaic_hankel/src ; make clean
	cd sage_output/test ; make clean
	cd magma_output/test ; make clean
	cd vec_lzz_p_extra/test ; make clean
	cd lzz_pX_middle_product/test ; make clean
	cd lzz_pX_CRT/test ; make clean
	cd ZZ_CRT/test ; make clean
	cd mat_lzz_pX_extra/test ; make clean
	cd bivariate/test ; make clean
	cd cauchy_like_geometric_special/test ; make clean
	cd mosaic_hankel/test ; make clean
	cd ZZ_pX_CRT/test ; make clean
