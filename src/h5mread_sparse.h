#ifndef _H5MREAD_SPARSE_H_
#define _H5MREAD_SPARSE_H_

#include "ChunkIterator.h"
#include <Rdefines.h>

SEXP _h5mread_sparse(
	ChunkIterator *chunk_iter,
	const size_t *ans_dim
);

#endif  /* _H5MREAD_SPARSE_H_ */

