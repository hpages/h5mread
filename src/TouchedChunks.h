#ifndef _TOUCHEDCHUNKS_H_
#define _TOUCHEDCHUNKS_H_

#include "H5DSetDescriptor.h"
#include "h5mread_helpers.h"
#include <Rdefines.h>
#include "S4Vectors_interface.h"

/* A data structure for representing the set of all the HDF5 chunks touched
   by user-supplied N-index. */
typedef struct all_tchunks_t {
	const H5DSetDescriptor *h5dset;
	SEXP index;                 /* user supplied N-index */
	LLongAEAE *breakpoint_bufs;
	LLongAEAE *tchunkidx_bufs;  /* 0-based indices of touched chunks
				       along each dim */
	size_t *num_tchunks;        /* nb of touched chunks along each dim */
	long long int total_num_tchunks;
} AllTChunks;

typedef struct tchunk_viewports_t {
	H5Viewport h5chunk_vp, mem_vp;
} TChunkViewports;

/* A data structure for iterating over the chunks of an HDF5 dataset. */
typedef struct tchunk_iterator_t {
	const AllTChunks *all_tchunks;
	size_t *tchunk_midx_buf;
	int moved_along;
	long long int tchunk_rank;
	TChunkViewports tchunk_vps;
} TChunkIterator;

/* A data structure for storing the data of a full chunk. */
typedef struct chunk_data_buffer_t {
	size_t data_length;
	hid_t data_type_id;
	size_t data_type_size;
	size_t data_size;
	hid_t data_space_id;
	void *data;
	H5Viewport data_vp;
	void *compressed_data;  /* experimental! */
} ChunkDataBuffer;

int _init_AllTChunks(
	AllTChunks *all_tchunks,
	const H5DSetDescriptor *h5dset,
	SEXP index,
	size_t *selection_dim
);

int _alloc_TChunkViewports(
	TChunkViewports *tchunk_vps,
	int ndim,
	int alloc_full_mem_vp
);

void _free_TChunkViewports(
	TChunkViewports *tchunk_vps
);

int _get_tchunk(
	const AllTChunks *all_tchunks,
	long long int i,
	size_t *tchunk_midx_buf,
	TChunkViewports *tchunk_vps
);

int _tchunk_is_truncated(
	const H5DSetDescriptor *h5dset,
	const H5Viewport *h5chunk_vp
);

int _tchunk_is_fully_selected(
	int ndim,
	const TChunkViewports *tchunk_vps
);

void _destroy_TChunkIterator(
	TChunkIterator *tchunk_iter
);

int _init_TChunkIterator(
	TChunkIterator *tchunk_iter,
	const AllTChunks *all_tchunks,
	int alloc_full_mem_vp
);

int _next_tchunk(
	TChunkIterator *tchunk_iter
);

void _print_tchunk_info(
	const TChunkIterator *tchunk_iter
);

void _destroy_ChunkDataBuffer(
	ChunkDataBuffer *chunk_data_buf
);

int _init_ChunkDataBuffer(
	ChunkDataBuffer *chunk_data_buf,
	const H5DSetDescriptor *h5dset,
	int use_Rtype
);

int _load_chunk(
	const H5DSetDescriptor *h5dset,
	const TChunkViewports *tchunk_vps,
	ChunkDataBuffer *chunk_data_buf,
	int use_H5Dread_chunk
);

int _reclaim_vlen_bufs(
	ChunkDataBuffer *chunk_data_buf
);

#endif  /* _TOUCHEDCHUNKS_H_ */

