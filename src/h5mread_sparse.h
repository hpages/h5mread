#ifndef _H5MREAD_SPARSE_H_
#define _H5MREAD_SPARSE_H_

#include "TouchedChunks.h"
#include <Rdefines.h>

SEXP _h5mread_sparse(
	const TouchedChunks *touched_chunks,
	const size_t *ans_dim
);

#endif  /* _H5MREAD_SPARSE_H_ */

