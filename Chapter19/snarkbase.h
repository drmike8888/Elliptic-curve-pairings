void p_i(mpz_t p, int i, mpz_t *list, int n);
void li_lj(mpz_t *coef, int i, int j, mpz_t *list, int n);
void liofx(mpz_t *coef, int i, mpz_t *list, int n);
void liljofx(mpz_t *coef, int i, int j, mpz_t *list, int n);
mpz_t* all_lilj(mpz_t *list, int n);
void lcalc(mpz_t rslt, mpz_t z, mpz_t *coef, int deg);
void tofzgrth(mpz_t t, mpz_t z, mpz_t *list, int n);
void matflat(mpz_t *vector, mpz_t *mat, int width, mpz_t *coef, int length);
