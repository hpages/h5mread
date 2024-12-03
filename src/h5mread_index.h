#ifndef _H5MREAD_INDEX_H_
#define _H5MREAD_INDEX_H_

#include "TouchedChunks.h"
#include <Rdefines.h>

SEXP _h5mread_index(
	const AllTChunks *all_tchunks,
	int method,
	int use_H5Dread_chunk,
	const size_t *ans_dim
);

#endif  /* _H5MREAD_INDEX_H_ */

