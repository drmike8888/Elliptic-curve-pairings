modulo.o:	modulo.c modulo.h
	gcc -c -o modulo.o modulo.c

poly.o:  	poly.h poly.c
	gcc -c -o poly.o poly.c

eliptic.o:	eliptic.c eliptic.h
	gcc -c -o eliptic.o eliptic.c

poly_eliptic.o: poly_eliptic.h poly_eliptic.c
	gcc -c -o poly_eliptic.o poly_eliptic.c

test_mod:	test_mod.c modulo.o poly.o
	gcc -o test_mod test_mod.c modulo.o poly.o -lgmp

test_irrd:	test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o
	gcc -o test_irrd test_irrd.c modulo.o poly.o eliptic.o poly_eliptic.o -lgmp

pairing_sweep_alpha:	pairing_sweep_alpha.c 
	gcc -o pairing_sweep_alpha pairing_sweep_alpha.c -lgmp -lm

get_curve:	get_curve.c modulo.o poly.o eliptic.o
	gcc -o get_curve get_curve.c modulo.o poly.o eliptic.o -lgmp

eliptic_test:	eliptic_test.c modulo.o eliptic.o
	gcc -o eliptic_test eliptic_test.c eliptic.o modulo.o -lgmp

pairing.o:	pairing.c pairing.h
	gcc -c -o pairing.o pairing.c

weil_6_bit_pairing:	weil_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o weil_6_bit_pairing weil_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

tate_6_bit_pairing:	tate_6_bit_pairing.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o tate_6_bit_pairing tate_6_bit_pairing.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

quotient_group:	quotient_group.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o quotient_group quotient_group.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

signature.o: signature.c signature.h
	gcc -c -o signature.o signature.c

signatures_11: signatures_11.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o signatures_11 signatures_11.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

signatures_11_keygen: signatures_11_keygen.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o signatures_11_keygen signatures_11_keygen.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

base_curves_embed:  base_curves_embed.c modulo.o
	gcc -o base_curves_embed base_curves_embed.c modulo.o -lgmp

base_protocols.o: base_protocols.c base_protocols.h
	gcc -c -o base_protocols.o base_protocols.c

base_test:  base_test.c base_protocols.o modulo.o eliptic.o 
	gcc -o base_test base_test.c base_protocols.o eliptic.o modulo.o -lgmp -lk12

base_curve_gen: base_curve_gen.c modulo.o eliptic.o
	gcc -o base_curve_gen base_curve_gen.c eliptic.o modulo.o -lgmp

poly_exp_test:	poly_exp_test.c poly.o modulo.o
	gcc -o poly_exp_test poly_exp_test.c poly.o modulo.o -lgmp

tiny_pwer_chk:	tiny_pwer_chk.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o
	gcc -o tiny_pwer_chk tiny_pwer_chk.c modulo.o poly.o eliptic.o poly_eliptic.o pairing.o -lgmp

snarkbase.o: snarkbase.c modulo.h snarkbase.h
	gcc -c -o snarkbase.o snarkbase.c

snark_qap: snark_qap.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o snark_qap snark_qap.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

snark_crs: snark_crs.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o snark_crs snark_crs.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

snark_proof: snark_proof.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o snark_proof snark_proof.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

snark_verify: snark_verify.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o snark_verify snark_verify.c snarkbase.o modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12

signatures_2_keygen: signatures_2_keygen.c modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o
	gcc -o signatures_2_keygen signatures_2_keygen.c  modulo.o eliptic.o poly.o poly_eliptic.o pairing.o signature.o -lgmp -lk12
