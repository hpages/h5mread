#ifndef _UASELECTION_H_
#define _UASELECTION_H_

#include <Rdefines.h>
#include "S4Vectors_interface.h"

/* Terminology:
   - uaselection: user-supplied array selection.
   - chips: the "chips" in the uaselection are its connected components i.e.
            its contiguous block-like components.
 */

/* Like VECTOR_ELT(x, i) except that 'x' can be R_NilValue. */
#define GET_LIST_ELT(x, i) ((x) != R_NilValue ? VECTOR_ELT(x, i) : R_NilValue)

int _shallow_check_uaselection(
	int ndim,
	SEXP starts,
	SEXP counts
);

long long int _check_uaselection(
	int ndim,
	const size_t *dim,
	SEXP starts,
	SEXP counts,
	size_t *uaselection_dim_buf
);

SEXP C_check_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

long long int _check_ordered_uaselection(
	int ndim,
	const size_t *dim,
	SEXP starts,
	SEXP counts,
	size_t *uaselection_dim_buf,
	size_t *nstart_buf,
	size_t *nchip_buf,
	size_t *last_chip_start_buf
);

SEXP C_check_ordered_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _uaselection_can_be_reduced(
	int ndim,
	SEXP starts,
	const size_t *nstart,
	const size_t *nchip
);

SEXP _reduce_uaselection(
	int ndim,
	SEXP starts, SEXP counts,
	const size_t *uaselection_dim,
	const size_t *nchip,
	const size_t *last_chip_start
);

SEXP C_reduce_uaselection(
	SEXP dim,
	SEXP starts,
	SEXP counts
);

int _map_starts_to_chunks(
	int ndim,
	const size_t *dim,
	const size_t *chunkdim,
	SEXP starts,
	size_t *nstart_buf,
	LLongAEAE *breakpoint_bufs,
	LLongAEAE *tchunkidx_bufs
);

SEXP C_map_starts_to_chunks(
	SEXP starts,
	SEXP dim,
	SEXP chunkdim
);


/****************************************************************************
 * Inline functions
 */

static inline long long int get_trusted_elt(SEXP x, R_xlen_t i)
{
	return IS_INTEGER(x) ? (long long int) INTEGER(x)[i] :
			       (long long int) REAL(x)[i];
}

static inline size_t *R_alloc0_size_t_array(int n)
{
        size_t *x = (size_t *) R_alloc((size_t) n, sizeof(size_t));
        memset(x, 0, sizeof(size_t) * n);
        return x;
}

#endif  /* _UASELECTION_H_ */

