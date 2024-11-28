/****************************************************************************
 *             Manipulation of a user-supplied array selection              *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "uaselection.h"

#include "global_errmsg_buf.h"

#include <limits.h>  /* for INT_MAX, LLONG_MAX, LLONG_MIN */
//#include <time.h>


/****************************************************************************
 * Low-level helpers
 */

static int check_INTEGER_or_NUMERIC(SEXP x, const char *what, int along)
{
	if (!(IS_INTEGER(x) || IS_NUMERIC(x))) {
		PRINT_TO_ERRMSG_BUF("'%s[[%d]]' must be an "
				    "integer vector (or NULL)",
				    what, along + 1);
		return -1;
	}
	return 0;
}

static int shallow_check_count(SEXP count, R_xlen_t n, int along)
{
	if (count == R_NilValue)
		return 0;
	if (check_INTEGER_or_NUMERIC(count, "counts", along) < 0)
		return -1;
	if (XLENGTH(count) != n) {
		PRINT_TO_ERRMSG_BUF("'starts[[%d]]' and 'counts[[%d]]' "
				    "must have the same length",
				    along + 1, along + 1);
		return -1;
	}
	return 0;
}

#define	NOT_A_FINITE_NUMBER(x) \
	(R_IsNA(x) || R_IsNaN(x) || (x) == R_PosInf || (x) == R_NegInf)

static inline int get_untrusted_elt(SEXP x, R_xlen_t i, long long int *val,
				    const char *what, int along)
{
	if (IS_INTEGER(x)) {
		int tmp1 = INTEGER(x)[i];
		if (tmp1 == NA_INTEGER) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%ld] is NA", what, i + 1);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%ld] is NA",
					    what, along + 1, i + 1);
		    return -1;
		}
		*val = (long long int) tmp1;
	} else {
		double tmp2 = REAL(x)[i];
		if (NOT_A_FINITE_NUMBER(tmp2)) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%ld] is NA or NaN "
					    "or not a finite number",
					    what, i + 1);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%ld] is NA or NaN "
					    "or not a finite number",
					    what, along + 1, i + 1);
		    return -1;
		}
		if (tmp2 > (double) LLONG_MAX || tmp2 < (double) LLONG_MIN) {
		    if (along < 0)
			PRINT_TO_ERRMSG_BUF("%s[%ld] is too large (= %e)",
					    what, i + 1, tmp2);
		    else
			PRINT_TO_ERRMSG_BUF("%s[[%d]][%ld] is too large (= %e)",
					    what, along + 1, i + 1, tmp2);
		    return -1;
		}
		*val = (long long int) tmp2;
	}
	return 0;
}

static inline void set_trusted_elt(SEXP x, R_xlen_t i, long long int val)
{
	if (IS_INTEGER(x)) {
		INTEGER(x)[i] = (int) val;
	} else {
		REAL(x)[i] = (double) val;
	}
	return;
}

/* Called at the very beginning of the various .Call entry points where
   it's used (and before any resource is allocated) so it's ok to error()
   immediately in case of error. */
static const size_t *check_dim(SEXP dim)
{
	if (!(IS_INTEGER(dim) || IS_NUMERIC(dim)))
		error("'dim' must be an integer vector");
	int ndim = LENGTH(dim);
	size_t *dim_p = R_alloc0_size_t_array(ndim);
	for (int along = 0; along < ndim; along++) {
		long long int d;
		int ret = get_untrusted_elt(dim, (R_xlen_t) along, &d,
					    "dim", -1);
		if (ret < 0)
			error("%s", _h5mread_global_errmsg_buf());
		if (d < 0)
			error("'dim' contains negative values");
		dim_p[along] = (size_t) d;
	}
	return dim_p;
}

static SEXP as_dim_vector(const size_t *uaselection_dim, int ndim)
{
	int use_doubles = 0;
	for (int along = 0; along < ndim; along++) {
		if (uaselection_dim[along] <= INT_MAX)
			continue;
		use_doubles = 1;
		break;
	}
	SEXP ans = PROTECT(allocVector(use_doubles ? REALSXP : INTSXP, ndim));
	for (int along = 0; along < ndim; along++) {
		size_t d = uaselection_dim[along];
		if (IS_INTEGER(ans)) {
			INTEGER(ans)[along] = (int) d;
		} else {
			REAL(ans)[along] = (double) d;
		}
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Shallow check of a user-supplied array selection
 */

/* Only check that each of 'starts' and 'counts' is either NULL or a list
   of length as 'ndim'.
   Return 0 if the uaselection is valid and -1 if it's not. */
int _shallow_check_uaselection(int ndim, SEXP starts, SEXP counts)
{
	if (starts == R_NilValue) {
		if (counts != R_NilValue) {
			PRINT_TO_ERRMSG_BUF(
				"'counts' must be NULL when 'starts' is NULL");
			return -1;
		}
		return 0;
	}
	if (!isVectorList(starts)) {  // IS_LIST() is broken
		PRINT_TO_ERRMSG_BUF("'starts' must be a list (or NULL)");
		return -1;
	}
	if (LENGTH(starts) != ndim) {
		PRINT_TO_ERRMSG_BUF(
			"Array has %d dimension%s but 'starts' has %d "
			"list element%s.\n  'starts' must have one "
			"list element per dimension in the dataset.",
			ndim, ndim > 1 ? "s" : "",
			LENGTH(starts), LENGTH(starts) > 1 ? "s" : "");
		return -1;
	}
	if (counts == R_NilValue)
		return 0;
	if (!isVectorList(counts)) {  // IS_LIST() is broken
		PRINT_TO_ERRMSG_BUF("'counts' must be a list (or NULL)");
		return -1;
	}
	if (LENGTH(counts) != ndim) {
		PRINT_TO_ERRMSG_BUF("'counts' must have one list element "
				    "per list element in 'starts'");
		return -1;
	}
	return 0;
}


/****************************************************************************
 * Deep check of a user-supplied array selection
 */

static void set_errmsg_for_uaselection_beyond_dim(
		int along1, R_xlen_t i,
		int no_counts)
{
	const char *msg = "selection must be within extent of "
			  "array, but you\n  have:";
	if (no_counts)
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%ld] "
			"> dimension %d in array",
			msg, along1, i + 1, along1);
	else
		PRINT_TO_ERRMSG_BUF(
			"%s starts[[%d]][%ld] + counts[[%d]][%ld] - 1 "
			"> dimension %d in array",
			msg, along1, i + 1, along1, i + 1, along1);
	return;
}

static void set_errmsg_for_non_strictly_ascending_uaselection(
		int along1, R_xlen_t i,
		int no_counts)
{
	const char *msg = "selection must be strictly ascending "
			  "along each dimension, but\n  you have:";
	if (no_counts)
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%ld] <= starts[[%d]][%ld]",
				    msg, along1, i + 1, along1, i);
	else
		PRINT_TO_ERRMSG_BUF("%s starts[[%d]][%ld] < starts[[%d]][%ld] "
				    " + counts[[%d]][%ld]",
				    msg, along1, i + 1, along1, i, along1, i);
	return;
}

static inline int get_untrusted_start(SEXP start, R_xlen_t i, long long int *s,
				      long long int min_start,
				      int along, int no_counts)
{
	if (get_untrusted_elt(start, i, s, "starts", along) < 0)
		return -1;
	if (*s < 1) {
		PRINT_TO_ERRMSG_BUF("starts[[%d]][%ld] is < 1",
				    along + 1, i + 1);
		return -1;
	}
	if (*s < min_start) {
		set_errmsg_for_non_strictly_ascending_uaselection(
			along + 1, i, no_counts);
		return -1;
	}
	return 0;
}

static long long int check_uaselection_along(int along,
		SEXP start, SEXP count, size_t d)
{
	if (start == R_NilValue) {
		if (count == R_NilValue)
			return (long long int) d;
		PRINT_TO_ERRMSG_BUF(
			"if 'starts[[%d]]' is NULL then 'counts' "
			"or 'counts[[%d]]' must also be NULL",
			along + 1, along + 1);
		return -1;
	}
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	R_xlen_t n = XLENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	/* Walk on the 'start' elements. */
	for (R_xlen_t i = 0; i < n; i++) {
		/* Last arg ('no_counts') is ignored when 4th arg ('min_start')
		   is set to 0. */
		long long int s;
		int ret = get_untrusted_start(start, i, &s, 0, along, 0);
		if (ret < 0)
			return -1;
		if (s > (long long int) d) {
			set_errmsg_for_uaselection_beyond_dim(
				along + 1, i, 1);
			return -1;
		}
	}
	if (count == R_NilValue)
		return (long long int) n;
	/* Walk on the 'count' (and 'start') elements. */
	long long int uaselection_dim = 0;
	for (R_xlen_t i = 0; i < n; i++) {
		long long int c;
		int ret = get_untrusted_elt(count, i, &c, "counts", along);
		if (ret < 0)
			return -1;
		if (c == 0)
			continue;
		if (c < 0) {
			PRINT_TO_ERRMSG_BUF("counts[[%d]][%ld] is < 0",
					    along + 1, i + 1);
			return -1;
		}
		long long int s = get_trusted_elt(start, i);
		long long int e = s + c - 1;
		if (e > (long long int) d) {
			set_errmsg_for_uaselection_beyond_dim(
				along + 1, i, 0);
			return -1;
		}
		uaselection_dim += c;
	}
	return uaselection_dim;
}

/* 'dim' must be an array of 'ndim' elements.

   'starts' and 'counts' are **assumed** to be NULL or lists of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   'uaselection_dim_buf' must point to an array of 'ndim' elements.

   To perform in-place replacement of 'dim' with 'uaselection_dim_buf'
   call _check_uaselection(ndim, dim_buf, starts, counts, dim_buf)
   where 'dim_buf' is a writable array containing the dimensions of the
   original array. */
long long int _check_uaselection(int ndim, const size_t *dim,
		SEXP starts, SEXP counts, size_t *uaselection_dim_buf)
{
	long long int uaselection_len = 1;
	for (int along = 0; along < ndim; along++) {
		SEXP start = GET_LIST_ELT(starts, along);
		SEXP count = GET_LIST_ELT(counts, along);
		long long int uaselection_dim =
				check_uaselection_along(along, start, count,
							dim[along]);
		if (uaselection_dim < 0)
			return -1;
		uaselection_dim_buf[along] = (size_t) uaselection_dim;
		uaselection_len *= uaselection_dim;
	}
	return uaselection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the user-supplied array selection.
 */
SEXP C_check_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const size_t *dim_p = check_dim(dim);
	int ndim = LENGTH(dim);
	int ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	size_t *uaselection_dim_buf = R_alloc0_size_t_array(ndim);
	long long int uaselection_len =
			_check_uaselection(ndim, dim_p, starts, counts,
					   uaselection_dim_buf);
	if (uaselection_len < 0)
		error("%s", _h5mread_global_errmsg_buf());
	return as_dim_vector(uaselection_dim_buf, ndim);
}


/****************************************************************************
 * Deep check of an ordered user-supplied array selection (in preparation for
 * its reduction)
 *
 * The "chips" in the user-supplied array selection are its connected
 * components i.e. its contiguous block-like components.
 */

static long long int check_ordered_uaselection_along_NULL_start(int along,
			SEXP count, size_t d,
			size_t *nstart_buf, size_t *nchip_buf,
			size_t *last_chip_start_buf)
{
	if (count != R_NilValue) {
		PRINT_TO_ERRMSG_BUF(
			"if 'starts[[%d]]' is NULL then 'counts' "
			"or 'counts[[%d]]' must also be NULL",
			along + 1, along + 1);
		return -1;
	}
	nstart_buf[along] = d;
	nchip_buf[along] = d != 0;
	last_chip_start_buf[along] = 1;
	return (long long int) d;
}

static long long int check_ordered_uaselection_along(int along,
			SEXP start, SEXP count, size_t d,
			size_t *nstart_buf, size_t *nchip_buf,
			size_t *last_chip_start_buf)
{
	if (start == R_NilValue)
		return check_ordered_uaselection_along_NULL_start(along,
				count, d,
				nstart_buf, nchip_buf, last_chip_start_buf);
	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;
	R_xlen_t n = XLENGTH(start);
	if (shallow_check_count(count, n, along) < 0)
		return -1;
	nstart_buf[along] = (size_t) n;
	nchip_buf[along] = 0;
	long long int uaselection_dim, min_start = 0;
	if (count == R_NilValue) {
		/* Walk on the 'start' elements. */
		for (R_xlen_t i = 0; i < n; i++) {
			long long int s;
			int ret = get_untrusted_start(start, i, &s, min_start,
						      along, 1);
			if (ret < 0)
				return -1;
			if (s != min_start) {
				nchip_buf[along]++;
				last_chip_start_buf[along] = s;
			}
			min_start = s + 1;
			if (s > (long long int) d) {
				set_errmsg_for_uaselection_beyond_dim(
					along + 1, i, 1);
				return -1;
			}
		}
		uaselection_dim = (long long int) n;
	} else {
		/* Walk on the 'start' and 'count' elements. */
		uaselection_dim = 0;
		for (R_xlen_t i = 0; i < n; i++) {
			long long int c;
			int ret = get_untrusted_elt(count, i, &c,
						    "counts", along);
			if (ret < 0)
				return -1;
			if (c < 0) {
				PRINT_TO_ERRMSG_BUF("counts[[%d]][%ld] is < 0",
						    along + 1, i + 1);
				return -1;
			}
			if (c == 0)
				continue;
			long long int s;
			ret = get_untrusted_start(start, i, &s, min_start,
						  along, 0);
			if (ret < 0)
				return -1;
			if (s != min_start) {
				nchip_buf[along]++;
				last_chip_start_buf[along] = s;
			}
			min_start = s + c;
			if (min_start - 1 > (long long int) d) {
				set_errmsg_for_uaselection_beyond_dim(
					along + 1, i, 0);
				return -1;
			}
			uaselection_dim += c;
		}
	}
	return uaselection_dim;
}

/* 'dim' must be an array of 'ndim' elements.

   'starts' and 'counts' are **assumed** to be NULL or lists of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   Each of 'uaselection_dim_buf', 'nstart_buf', 'nchip_buf', and
   'last_chip_start_buf' must point to an array of 'ndim' elements.

   To perform in-place replacement of 'dim' with 'uaselection_dim_buf'
   call _check_ordered_uaselection(ndim, dim_buf, starts, counts, dim_buf, ...)
   where 'dim_buf' is a writable array containing the dimensions of the
   original array. */
long long int _check_ordered_uaselection(int ndim, const size_t *dim,
			SEXP starts, SEXP counts,
			size_t *uaselection_dim_buf,
			size_t *nstart_buf, size_t *nchip_buf,
			size_t *last_chip_start_buf)
{
	long long int uaselection_len = 1;
	for (int along = 0; along < ndim; along++) {
		SEXP start = GET_LIST_ELT(starts, along);
		SEXP count = GET_LIST_ELT(counts, along);
		long long int uaselection_dim =
			check_ordered_uaselection_along(along,
					start, count,
					dim[along],
					nstart_buf, nchip_buf,
					last_chip_start_buf);
		if (uaselection_dim < 0)
			return -1;
		uaselection_dim_buf[along] = (size_t) uaselection_dim;
		uaselection_len *= uaselection_dim;
	}
	return uaselection_len;
}

/* --- .Call ENTRY POINT ---
 * Return the dimensions of the user-supplied array selection.
 */
SEXP C_check_ordered_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const size_t *dim_p = check_dim(dim);
	int ndim = LENGTH(dim);
	int ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	size_t *uaselection_dim_buf = R_alloc0_size_t_array(ndim);
	size_t *nstart_buf          = R_alloc0_size_t_array(ndim);
	size_t *nchip_buf           = R_alloc0_size_t_array(ndim);
	size_t *last_chip_start_buf = R_alloc0_size_t_array(ndim);
	long long int uaselection_len =
			_check_ordered_uaselection(ndim, dim_p,
					starts, counts,
					uaselection_dim_buf,
					nstart_buf, nchip_buf,
					last_chip_start_buf);
	if (uaselection_len < 0)
		error("%s", _h5mread_global_errmsg_buf());
	return as_dim_vector(uaselection_dim_buf, ndim);
}


/****************************************************************************
 * Reduce the user-supplied array selection
 */

int _uaselection_can_be_reduced(int ndim, SEXP starts,
		const size_t *nstart, const size_t *nchip)
{
	if (starts == R_NilValue)
		return 0;
	for (int along = 0; along < ndim; along++) {
		if (VECTOR_ELT(starts, along) == R_NilValue)
			continue;
		/* nchip[along] should always be <= nstart[along] */
		if (nchip[along] < nstart[along])
			return 1;
	}
	return 0;
}

static SEXP dup_or_coerce_to_INTSXP(SEXP x, int dup)
{
	if (dup)
		return duplicate(x);
	R_xlen_t x_len = XLENGTH(x);
	SEXP ans = PROTECT(NEW_INTEGER(x_len));
	for (R_xlen_t i = 0; i < x_len; i++)
		INTEGER(ans)[i] = (int) REAL(x)[i];
	UNPROTECT(1);
	return ans;
}

/*
 * Note that this does something similar to what coercion from integer (or
 * numeric) to IRanges does (see .Call entry point "IRanges_from_integer"
 * in IRanges). However we cannot re-use this here because we want to be able
 * to handle start values that are >= 2^31 which this coercion doesn't support
 * at the moment.
 */
static void stitch_uaselection(SEXP start_in, SEXP count_in,
			       SEXP start_out, int *count_out)
{
	R_xlen_t n = XLENGTH(start_in);
	long long int min_start = 0;
	R_xlen_t j = -1;
	if (count_in == R_NilValue) {
		for (R_xlen_t i = 0; i < n; i++) {
			long long int s = get_trusted_elt(start_in, i);
			if (s != min_start) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = 1;
			} else {
				/* FIXME: Risk of integer overflow! */
				count_out[j]++;
			}
			min_start = s + 1;
		}
	} else {
		for (R_xlen_t i = 0; i < n; i++) {
			long long int c = get_trusted_elt(count_in, i);
			if (c == 0)
				continue;
			long long int s = get_trusted_elt(start_in, i);
			if (s != min_start) {
				j++;
				set_trusted_elt(start_out, j, s);
				count_out[j] = c;
			} else {
				/* FIXME: Risk of integer overflow! */
				count_out[j] += c;
			}
			min_start = s + c;
		}
	}
	return;
}

static void reduce_uaselection_along(int along,
				    SEXP start, SEXP count,
				    const size_t *uaselection_dim,
				    const size_t *nchip,
				    const size_t *last_chip_start,
				    SEXP reduced_starts, SEXP reduced_counts)
{
	size_t n = XLENGTH(start);
	if (nchip[along] == n) {
		/* Nothing to stitch. */
		int dup = IS_INTEGER(start) || last_chip_start[along] > INT_MAX;
		SEXP reduced_start =
			PROTECT(dup_or_coerce_to_INTSXP(start, dup));
		SET_VECTOR_ELT(reduced_starts, along, reduced_start);
		UNPROTECT(1);
		if (uaselection_dim[along] == n)
			return;
		dup = IS_INTEGER(count);
		SEXP reduced_count =
			PROTECT(dup_or_coerce_to_INTSXP(count, dup));
		SET_VECTOR_ELT(reduced_counts, along, reduced_count);
		UNPROTECT(1);
		return;
	}
	/* Stitch. */
	SEXPTYPE type = last_chip_start[along] <= INT_MAX ? INTSXP : REALSXP;
	SEXP reduced_start =
		PROTECT(allocVector(type, (R_xlen_t) nchip[along]));
	SET_VECTOR_ELT(reduced_starts, along, reduced_start);
	UNPROTECT(1);
	SEXP reduced_count =
		PROTECT(NEW_INTEGER((R_xlen_t) nchip[along]));
	SET_VECTOR_ELT(reduced_counts, along, reduced_count);
	UNPROTECT(1);
	stitch_uaselection(start, count, reduced_start, INTEGER(reduced_count));
	return;
}

SEXP _reduce_uaselection(int ndim, SEXP starts, SEXP counts,
			 const size_t *uaselection_dim,
			 const size_t *nchip,
			 const size_t *last_chip_start)
{
	//clock_t t0 = clock();
	SEXP reduced_starts = PROTECT(NEW_LIST(ndim));
	SEXP reduced_counts = PROTECT(NEW_LIST(ndim));
	if (starts != R_NilValue) {
		for (int along = 0; along < ndim; along++) {
			SEXP start = VECTOR_ELT(starts, along);
			if (start == R_NilValue)
				continue;
			SEXP count = GET_LIST_ELT(counts, along);
			reduce_uaselection_along(along,
					start, count,
					uaselection_dim,
					nchip, last_chip_start,
					reduced_starts, reduced_counts);
		}
	}
	SEXP ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, reduced_starts);
	SET_VECTOR_ELT(ans, 1, reduced_counts);
	UNPROTECT(3);
	//printf("time 2nd pass: %e\n", (1.0 * clock() - t0) / CLOCKS_PER_SEC);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Negative values in 'dim' are treated as infinite dimensions.
 * Return a list of length 2 or NULL if the user-supplied array selection
 * could not be reduced.
 * When returning a list of length 2:
 *   - The 1st list element is the list of reduced starts.
 *   - The 2nd list element is the list of reduced counts.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_reduce_uaselection(SEXP dim, SEXP starts, SEXP counts)
{
	const size_t *dim_p = check_dim(dim);
	int ndim = LENGTH(dim);
	int ret = _shallow_check_uaselection(ndim, starts, counts);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	size_t *uaselection_dim_buf = R_alloc0_size_t_array(ndim);
	size_t *nstart_buf          = R_alloc0_size_t_array(ndim);
	size_t *nchip_buf           = R_alloc0_size_t_array(ndim);
	size_t *last_chip_start_buf = R_alloc0_size_t_array(ndim);

	/* 1st pass */
	long long int uaselection_len =
			_check_ordered_uaselection(ndim, dim_p,
					starts, counts,
					uaselection_dim_buf,
					nstart_buf, nchip_buf,
					last_chip_start_buf);
	if (uaselection_len < 0)
		error("%s", _h5mread_global_errmsg_buf());
	if (!_uaselection_can_be_reduced(ndim, starts, nstart_buf, nchip_buf))
		return R_NilValue;

	/* 2nd pass */
	return _reduce_uaselection(ndim, starts, counts,
				   uaselection_dim_buf, nchip_buf,
				   last_chip_start_buf);
}


/****************************************************************************
 * Map the user-supplied array selection to the physical chunks
 */

static int map_start_to_chunks(int along,
		size_t d, size_t chunkd, SEXP start,
		size_t *nstart_buf,
		LLongAE *breakpoint_buf, LLongAE *tchunkidx_buf)
{
	if (start == R_NilValue) {
		if (nstart_buf != NULL)
			nstart_buf[along] = d;
		return 0;
	}

	if (check_INTEGER_or_NUMERIC(start, "starts", along) < 0)
		return -1;

	if (LLongAE_get_nelt(breakpoint_buf) != 0 ||
	    LLongAE_get_nelt(tchunkidx_buf) != 0)
	{
		/* Should never happen! */
		PRINT_TO_ERRMSG_BUF("internal error: map_start_to_chunks() "
				    "was called with non-empty breakpoint "
				    "or tchunkidx buffers");
		return -1;
	}

	R_xlen_t n = XLENGTH(start);
	if (nstart_buf != NULL)
		nstart_buf[along] = n;

	if (n == 0)
		return 0;

	/* Get 's' and 'tchunkidx' for 1st 'start' element. */
	long long int s;
	int ret = get_untrusted_start(start, 0, &s, 1, along, 1);
	if (ret < 0)
		return -1;
	if (s > (long long int) d) {
		set_errmsg_for_uaselection_beyond_dim(along + 1, 0, 1);
		return -1;
	}
	long long int tchunkidx = (s - 1) / chunkd;

	/* Walk on the remaining 'start' elements. */
	size_t ntchunk = 0;
	for (R_xlen_t i = 1; i < n; i++) {
		long long int min_start = s + 1;
		ret = get_untrusted_start(start, i, &s, min_start, along, 1);
		if (ret < 0)
			return -1;
		if (s > (long long int) d) {
			set_errmsg_for_uaselection_beyond_dim(along + 1, i, 1);
			return -1;
		}
		long long int prev_tchunkidx = tchunkidx;
		tchunkidx = (s - 1) / chunkd;
		if (tchunkidx > prev_tchunkidx) {
			LLongAE_insert_at(breakpoint_buf, ntchunk,
					  (long long int) i);
			LLongAE_insert_at(tchunkidx_buf, ntchunk,
					  prev_tchunkidx);
			ntchunk++;
		}
	}
	LLongAE_insert_at(breakpoint_buf, ntchunk, (long long int) n);
	LLongAE_insert_at(tchunkidx_buf, ntchunk, tchunkidx);
	return 0;
}

/* 'dim', 'chunkdim', and 'nstart_buf' must point to arrays of 'ndim'
   elements (alternatively 'nstart_buf' can be set to NULL).

   'starts' is **assumed** to be NULL or a list of length 'ndim'.
   This should have been already checked by _shallow_check_uaselection() so is
   not checked again.

   'breakpoint_bufs' and 'tchunkidx_bufs' must be of length 'ndim'. */
int _map_starts_to_chunks(int ndim, const size_t *dim, const size_t *chunkdim,
		SEXP starts, size_t *nstart_buf,
		LLongAEAE *breakpoint_bufs, LLongAEAE *tchunkidx_bufs)
{
	for (int along = 0; along < ndim; along++) {
		SEXP start = GET_LIST_ELT(starts, along);
		int ret = map_start_to_chunks(along,
					dim[along], chunkdim[along],
					start, nstart_buf,
					breakpoint_bufs->elts[along],
					tchunkidx_bufs->elts[along]);
		if (ret < 0)
			return -1;
	}
	return 0;
}

static SEXP to_numeric_LIST(int ndim, const LLongAEAE *aeae, SEXP starts)
{
	SEXP ans = PROTECT(NEW_LIST(ndim));
	if (starts != R_NilValue) {
		for (int along = 0; along < ndim; along++) {
			SEXP start = VECTOR_ELT(starts, along);
			if (start == R_NilValue)
				continue;
			const LLongAE *ae = aeae->elts[along];
			R_xlen_t ans_elt_len = LLongAE_get_nelt(ae);
			SEXP ans_elt = PROTECT(NEW_NUMERIC(ans_elt_len));
			for (R_xlen_t i = 0; i < ans_elt_len; i++)
				REAL(ans_elt)[i] = (double) ae->elts[i];
			SET_VECTOR_ELT(ans, along, ans_elt);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT ---
 * Return a list of length 2:
 *   - The 1st list element is the list of break points along each dim.
 *   - The 2nd list element is the list of touched chunk ids along each dim.
 * The 2 lists have the same length as 'starts'. Also they have the same
 * shape (i.e. same lengths()).
 */
SEXP C_map_starts_to_chunks(SEXP starts, SEXP dim, SEXP chunkdim)
{
	const size_t *dim_p = check_dim(dim);
	int ndim = LENGTH(dim);
	int ret = _shallow_check_uaselection(ndim, starts, R_NilValue);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	if (!(IS_INTEGER(chunkdim) || IS_NUMERIC(chunkdim)))
		error("'chunkdim' must be an integer vector");
	if (LENGTH(chunkdim) != ndim)
		error("'chunkdim' must have the same length as 'dim'");
	size_t *chunkdim_buf = R_alloc0_size_t_array(ndim);

	for (int i = 0; i < ndim; i++) {
		long long int chunkd;
		ret = get_untrusted_elt(chunkdim, (R_xlen_t) i, &chunkd,
					"chunkdim", -1);
		if (ret < 0)
			error("%s", _h5mread_global_errmsg_buf());
		if (chunkd < 0)
			error("'chunkdim' cannot contain negative values");
		if (chunkd == 0 && dim_p[i] != 0)
			error("values in 'chunkdim' cannot be 0 unless "
			      "their corresponding value\n  in 'dim' is "
			      "also 0");
		chunkdim_buf[i] = (size_t) chunkd;
	}

	size_t *nstart_buf = R_alloc0_size_t_array(ndim);
	LLongAEAE *breakpoint_bufs = new_LLongAEAE(ndim, ndim);
	LLongAEAE *tchunkidx_bufs  = new_LLongAEAE(ndim, ndim);
	ret = _map_starts_to_chunks(ndim, dim_p, chunkdim_buf,
			starts, nstart_buf,
			breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0)
		error("%s", _h5mread_global_errmsg_buf());

	SEXP ans = PROTECT(NEW_LIST(2));
	SEXP ans_elt = PROTECT(to_numeric_LIST(ndim, breakpoint_bufs, starts));
	SET_VECTOR_ELT(ans, 0, ans_elt);
	UNPROTECT(1);
	ans_elt = PROTECT(to_numeric_LIST(ndim, tchunkidx_bufs, starts));
	SET_VECTOR_ELT(ans, 1, ans_elt);
	UNPROTECT(2);
	return ans;
}

