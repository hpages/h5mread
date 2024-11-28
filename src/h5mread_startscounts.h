#ifndef _H5MREAD_STARTSCOUNTS_H_
#define _H5MREAD_STARTSCOUNTS_H_

#include "H5DSetDescriptor.h"
#include <Rdefines.h>

SEXP _compute_startscounts_ans_dim(
	const H5DSetDescriptor *h5dset,
	SEXP starts,
	SEXP counts,
	int noreduce,
	int method,
	size_t *ans_dim
);

SEXP _h5mread_startscounts(
	const H5DSetDescriptor *h5dset,
	SEXP startscounts,
	int noreduce,
	int method,
	const size_t *ans_dim
);

#endif  /* _H5MREAD_STARTSCOUNTS_H_ */

