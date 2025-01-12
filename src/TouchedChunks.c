/****************************************************************************
 *           Manipulating the "touched" chunks of an HDF5 dataset           *
 *                          and iterating over them                         *
 *                                 --------                                 *
 *                            Author: H. Pag\`es                            *
 ****************************************************************************/
#include "TouchedChunks.h"

#include "global_errmsg_buf.h"
#include "uaselection.h"

#include <stdlib.h>  /* for malloc, free */
#include <zlib.h>  /* for uncompress(), Z_OK, Z_MEM_ERROR, etc.. */


/****************************************************************************
 * Low-level helpers (non-exported)
 */

static int alloc_h5chunk_vp_mem_vp(int ndim,
		H5Viewport *h5chunk_vp,
		H5Viewport *mem_vp, int mem_vp_mode)
{
	if (_alloc_H5Viewport(h5chunk_vp, ndim, ALLOC_H5OFF_AND_H5DIM) < 0)
		return -1;
	if (_alloc_H5Viewport(mem_vp, ndim, mem_vp_mode) < 0) {
		_free_H5Viewport(h5chunk_vp);
		return -1;
	}
	return 0;
}

static void free_h5chunk_vp_mem_vp(H5Viewport *h5chunk_vp, H5Viewport *mem_vp)
{
	_free_H5Viewport(mem_vp);
	_free_H5Viewport(h5chunk_vp);
	return;
}

static int map_starts_to_h5chunks(const H5DSetDescriptor *h5dset,
		SEXP starts, size_t *nstart_buf,
		LLongAEAE *breakpoint_bufs, LLongAEAE *tchunkidx_bufs)
{
	int ndim = h5dset->ndim;
	size_t *dim_buf      = R_alloc0_size_t_array(ndim);
	size_t *chunkdim_buf = R_alloc0_size_t_array(ndim);
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		dim_buf[along]      = (size_t) h5dset->h5dim[h5along];
		chunkdim_buf[along] = (size_t) h5dset->h5chunkdim[h5along];
	}
	return _map_starts_to_chunks(ndim, dim_buf, chunkdim_buf,
				     starts, nstart_buf,
				     breakpoint_bufs, tchunkidx_bufs);
}

static long long int set_num_tchunks(const H5DSetDescriptor *h5dset,
		const SEXP starts,
		const LLongAEAE *tchunkidx_bufs,
		size_t *num_tchunks_buf)
{
	int ndim = h5dset->ndim;
	long long int total_num_tchunks = 1;  /* total nb of touched chunks */
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		SEXP start = GET_LIST_ELT(starts, along);
		size_t n;
		if (start != R_NilValue) {
			n = LLongAE_get_nelt(tchunkidx_bufs->elts[along]);
		} else {
			n = h5dset->h5nchunk[h5along];
		}
		total_num_tchunks *= num_tchunks_buf[along] = n;
	}
	return total_num_tchunks;
}

static void update_h5chunk_vp(const H5DSetDescriptor *h5dset,
		const size_t *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *tchunkidx_bufs,
		H5Viewport *h5chunk_vp)
{
	int ndim = h5dset->ndim;
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		size_t i = tchunk_midx[along];
		SEXP start = GET_LIST_ELT(starts, along);
		long long int tchunkidx;
		if (start != R_NilValue) {
			tchunkidx = tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = (long long int) i;
		}
		hsize_t chunkd = h5dset->h5chunkdim[h5along];
		hsize_t off = tchunkidx * chunkd;
		hsize_t d = h5dset->h5dim[h5along] - off;
		if (d > chunkd)
			d = chunkd;
		h5chunk_vp->h5off[h5along] = off;
		h5chunk_vp->h5dim[h5along] = d;
	}
	//Rprintf("# h5chunk_vp->h5off:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      Rprintf(" %llu", h5chunk_vp->h5off[h5along]);
	//Rprintf("\n");
	//Rprintf("# h5chunk_vp->h5dim:");
	//for (h5along = ndim - 1; h5along >= 0; h5along--)
	//      Rprintf(" %llu", h5chunk_vp->h5dim[h5along]);
	//Rprintf("\n");
	return;
}

static void update_mem_vp(const H5DSetDescriptor *h5dset,
		const size_t *tchunk_midx, int moved_along,
		SEXP starts, const LLongAEAE *breakpoint_bufs,
		const H5Viewport *h5chunk_vp, H5Viewport *mem_vp)
{
	int ndim = h5dset->ndim;
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		if (along > moved_along)
			break;
		size_t i = tchunk_midx[along];
		SEXP start = GET_LIST_ELT(starts, along);
		long long int off, d;
		if (start != R_NilValue ) {
			const long long int *breakpoint =
					breakpoint_bufs->elts[along]->elts;
			off = i == 0 ? 0 : breakpoint[i - 1];
			d = breakpoint[i] - off;
		} else {
			off = h5chunk_vp->h5off[h5along];
			d = h5chunk_vp->h5dim[h5along];
		}
		if (mem_vp->h5off != NULL) {
			mem_vp->h5off[h5along] = (hsize_t) off;
			mem_vp->h5dim[h5along] = (hsize_t) d;
		}
		mem_vp->off[along] = (size_t) off;
		mem_vp->dim[along] = (size_t) d;
	}
	//Rprintf("# mem_vp (offsets):");
	//for (along = 0; along < ndim; along++)
	//      Rprintf(" %lu", mem_vp->off[along]);
	//Rprintf("\n");
	//Rprintf("# mem_vp (dims):");
	//for (along = 0; along < ndim; along++)
	//      Rprintf(" %lu", mem_vp->dim[along]);
	//Rprintf("\n");
	return;
}

static void update_TChunkViewports(const H5DSetDescriptor *h5dset,
		const size_t *tchunk_midx, int moved_along,
		SEXP starts,
		const LLongAEAE *breakpoint_bufs,
		const LLongAEAE *tchunkidx_bufs,
		TChunkViewports *tchunk_vps)
{
	update_h5chunk_vp(h5dset,
			tchunk_midx, moved_along,
			starts, tchunkidx_bufs,
			&tchunk_vps->h5chunk_vp);
	update_mem_vp(h5dset,
			tchunk_midx, moved_along,
			starts, breakpoint_bufs,
			&tchunk_vps->h5chunk_vp,
			&tchunk_vps->mem_vp);
	return;
}


/****************************************************************************
 * read_h5chunk()
 *
 * Based on H5Dread_chunk(), which is NOT listed here for some mysterious
 * reasons: https://support.hdfgroup.org/HDF5/doc/RM/RM_H5D.html
 *
 * Header file for declaration: hdf5-1.10.3/src/H5Dpublic.h
 *
 * See hdf5-1.10.3/test/direct_chunk.c for plenty of examples.
 *
 * Call stack for H5Dread_chunk()
 *   H5Dread_chunk                (H5Dio.c)
 *     H5D__chunk_direct_read     (H5Dchunk.c)
 *       H5F_block_read           (H5Fio.c)
 *         H5PB_read              (H5PB.c)
 *           H5F__accum_read      (H5Faccum.c)
 *             or
 *           H5FD_read            (H5FDint.c)
 *             ??
 *
 * Call stack for H5Dread()
 *   H5Dread                      (H5Dio.c)
 *     H5D__read                  (H5Dio.c)
 *       H5D__chunk_read          (H5Dchunk.c)
 *         H5D__select_read
 *           or
 *         H5D__scatgath_read     (H5Dscatgath.c)
 *           H5D__gather_file     (H5Dscatgath.c)
 *       call ser_read member of a H5D_layout_ops_t object
 *            ??
 */

static int uncompress_chunk_data(const void *compressed_chunk_data,
				 hsize_t compressed_size,
				 void *uncompressed_chunk_data,
				 size_t uncompressed_size)
{
	uLong destLen = (uLong) uncompressed_size;
	int ret = uncompress((Bytef *) uncompressed_chunk_data, &destLen,
			 compressed_chunk_data, (uLong) compressed_size);
	if (ret == Z_OK) {
		if (destLen == uncompressed_size)
			return 0;
		PRINT_TO_ERRMSG_BUF("error in uncompress_chunk_data(): "
				    "chunk data smaller than expected "
				    "after decompression");
		return -1;
	}
	switch (ret) {
	    case Z_MEM_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough memory to uncompress chunk");
		break;
	    case Z_BUF_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "not enough room in output buffer");
		break;
	    case Z_DATA_ERROR:
		PRINT_TO_ERRMSG_BUF("error in uncompress(): "
				    "chunk data corrupted or incomplete");
		break;
	    default:
		PRINT_TO_ERRMSG_BUF("unknown error in uncompress()");
	}
	return -1;
}

static void transpose_bytes(const char *in, size_t nrow, size_t ncol, char *out)
{
	for (size_t i = 0; i < nrow; i++) {
		size_t in_offset = i;
		for (size_t j = 0; j < ncol; j++) {
			*(out++) = *(in + in_offset);
			in_offset += nrow;
		}
	}
	return;
}

/*
static void print_chunk_data(void *data, size_t data_length, size_t data_size)
{
	Rprintf("chunk data:");
	//for (size_t i = 0; i < data_size; i++) {
	//	if (i % 12 == 0)
	//		Rprintf("\n ");
	//	Rprintf(" '%c'", ((char *) data)[i]);
	//}
	for (size_t i = 0; i < data_length; i++) {
		if (i % 12 == 0)
			Rprintf("\n ");
		Rprintf(" %4d", ((int *) data)[i]);
	}
	Rprintf("\n");
	return;
}
*/

#define	CHUNK_COMPRESSION_OVERHEAD 8  // empirical (increase if necessary)

/* WARNING: read_h5chunk() is not ready yet! It is NOT working properly on
   some datasets:
      library(h5mread)
      library(ExperimentHub)
      hub <- ExperimentHub()
      fname0 <- hub[["EH1039"]]
      h5mread(fname0, "mm10/barcodes", list(1), method=4L)
      # [1] "AAACCTGAGATAGGAG-1"
      h5mread(fname0, "mm10/barcodes", list(1),
              method=4L, use.H5Dread_chunk=TRUE)
      # [1] "AAAAAAAAAAAAAAAAAAAA"
   Looks like the chunk data has been shuffled (transposed in that case)
   before being written to disk in order to improve compression.
   TODO: Investigate this further. I suspect we need to check whether a
   "Data shuffling filter" (H5Z_FILTER_SHUFFLE) was used at creation time.
   Check H5Pget_filter() for how to know whether this filter was used or not.
   There should be a way to retrieve information about how the data was
   shuffled. */
static int read_h5chunk(hid_t dset_id,
		const H5Viewport *h5chunk_vp,
		ChunkDataBuffer *chunk_data_buf)
{
	int ret;
	hsize_t chunk_storage_size;
	uint32_t filters;

	ret = H5Dget_chunk_storage_size(dset_id,
					h5chunk_vp->h5off,
					&chunk_storage_size);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dget_chunk_storage_size() "
				    "returned an error");
		return -1;
	}
	if (chunk_storage_size > chunk_data_buf->data_size +
				 CHUNK_COMPRESSION_OVERHEAD)
	{
		PRINT_TO_ERRMSG_BUF("chunk storage size (%llu) bigger "
				    "than expected (%llu + %d)",
				    chunk_storage_size,
				    (long long unsigned)
					chunk_data_buf->data_size,
				    CHUNK_COMPRESSION_OVERHEAD);
		return -1;
	}
	ret = H5Dread_chunk(dset_id, H5P_DEFAULT,
			    h5chunk_vp->h5off, &filters,
			    chunk_data_buf->compressed_data);
	if (ret < 0) {
		PRINT_TO_ERRMSG_BUF("H5Dread_chunk() returned an error");
		return -1;
	}

	//Rprintf("filters = %u\n", filters);

	//FIXME: This will error if chunk data is not compressed!
	//TODO: Decompress only if chunk data is compressed. There should be
	//a bit in the returned 'filters' that indicates this.
	ret = uncompress_chunk_data(chunk_data_buf->compressed_data,
				    chunk_storage_size,
				    chunk_data_buf->data,
				    chunk_data_buf->data_size);
	if (ret < 0)
		return -1;
	transpose_bytes(chunk_data_buf->data,
			chunk_data_buf->data_length,
			chunk_data_buf->data_type_size,
			chunk_data_buf->compressed_data);
	//print_chunk_data(chunk_data_buf->compressed_data,
	//		   chunk_data_buf->data_length,
	//		   chunk_data_buf->data_size);
	return 0;
}


/****************************************************************************
 * _init_AllTChunks()
 */

int _init_AllTChunks(AllTChunks *all_tchunks,
		const H5DSetDescriptor *h5dset, SEXP index,
		size_t *selection_dim)
{
	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	/* Set struct members 'h5dset' and 'index'. */
	all_tchunks->h5dset = h5dset;
	all_tchunks->index = index;

	/* Set struct members 'breakpoint_bufs' and 'tchunkidx_bufs'.
	   Also populate 'selection_dim' if not set to NULL. */
	int ndim = h5dset->ndim;
	LLongAEAE *breakpoint_bufs = new_LLongAEAE(ndim, ndim);
	LLongAEAE *tchunkidx_bufs  = new_LLongAEAE(ndim, ndim);
	int ret = map_starts_to_h5chunks(h5dset, index, selection_dim,
					 breakpoint_bufs, tchunkidx_bufs);
	if (ret < 0)
		return -1;
	all_tchunks->breakpoint_bufs = breakpoint_bufs;
	all_tchunks->tchunkidx_bufs  = tchunkidx_bufs;

	/* Set struct members 'num_tchunks' and 'total_num_tchunks'. */
	all_tchunks->num_tchunks = R_alloc0_size_t_array(ndim);
	all_tchunks->total_num_tchunks = set_num_tchunks(h5dset, index,
						tchunkidx_bufs,
						all_tchunks->num_tchunks);
	return 0;
}


/****************************************************************************
 * TChunkViewports utilities
 */

int _alloc_TChunkViewports(TChunkViewports *tchunk_vps,
			   int ndim, int alloc_full_mem_vp)
{
	return alloc_h5chunk_vp_mem_vp(ndim,
				       &tchunk_vps->h5chunk_vp,
				       &tchunk_vps->mem_vp,
				       alloc_full_mem_vp ? ALLOC_ALL_FIELDS
							 : ALLOC_OFF_AND_DIM);
}

void _free_TChunkViewports(TChunkViewports *tchunk_vps)
{
	free_h5chunk_vp_mem_vp(&tchunk_vps->h5chunk_vp, &tchunk_vps->mem_vp);
	return;
}

int _get_tchunk(const AllTChunks *all_tchunks, long long int i,
		size_t *tchunk_midx_buf,
		TChunkViewports *tchunk_vps)
{
	const H5DSetDescriptor *h5dset = all_tchunks->h5dset;
	//Rprintf("_get_tchunk(): Touched block %3lld/%3lld: ",
	//        i, all_tchunks->total_num_tchunks);
	for (int along = 0; along < h5dset->ndim; along++) {
		size_t n = all_tchunks->num_tchunks[along];
		tchunk_midx_buf[along] = i % n;
		//Rprintf(" %3ld/%ld", tchunk_midx_buf[along], n);
		i /= n;
	}
	//Rprintf("\n");
	if (i != 0) {
		PRINT_TO_ERRMSG_BUF("i >= total_num_tchunks");
		return -1;
	}
	update_TChunkViewports(h5dset,
			       tchunk_midx_buf,
			       h5dset->ndim,
			       all_tchunks->index,
			       all_tchunks->breakpoint_bufs,
			       all_tchunks->tchunkidx_bufs,
			       tchunk_vps);
	return 0;
}

/* Return 1 if the chunk that 'h5chunk_vp' is pointing at is "truncated"
   (a.k.a. "partial edge chunk" in HDF5's terminology), and 0 otherwise
   (i.e. if the new chunk is a full-size chunk). */
int _tchunk_is_truncated(const H5DSetDescriptor *h5dset,
			 const H5Viewport *h5chunk_vp)
{
	int ndim, h5along;
	hsize_t chunkd, d;

	ndim = h5dset->ndim;
	for (h5along = 0; h5along < ndim; h5along++) {
		chunkd = h5dset->h5chunkdim[h5along];
		d = h5chunk_vp->h5dim[h5along];
		if (d != chunkd)
			return 1;
	}
	return 0;
}

int _tchunk_is_fully_selected(int ndim, const TChunkViewports *tchunk_vps)
{
	int along, h5along, not_fully;

	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		not_fully = tchunk_vps->h5chunk_vp.h5dim[h5along] !=
			    (hsize_t) tchunk_vps->mem_vp.dim[along];
		if (not_fully)
			return 0;
	}
	return 1;
}


/****************************************************************************
 * TChunkIterator utilities
 */

void _destroy_TChunkIterator(TChunkIterator *tchunk_iter)
{
	_free_TChunkViewports(&tchunk_iter->tchunk_vps);
	return;
}

int _init_TChunkIterator(TChunkIterator *tchunk_iter,
			 const AllTChunks *all_tchunks,
			 int alloc_full_mem_vp)
{
	/* Set struct member 'all_tchunks'. */
	tchunk_iter->all_tchunks = all_tchunks;

	/* Allocate struct member 'tchunk_vps'. */
	int ndim = all_tchunks->h5dset->ndim;
	int ret = _alloc_TChunkViewports(&tchunk_iter->tchunk_vps,
					 ndim,
					 alloc_full_mem_vp);
	if (ret < 0)
		goto on_error;

	/* Set struct member 'tchunk_midx_buf'. */
	tchunk_iter->tchunk_midx_buf = R_alloc0_size_t_array(ndim);

	/* Set struct member 'tchunk_rank'. */
	tchunk_iter->tchunk_rank = -1;
	return 0;

    on_error:
	_destroy_TChunkIterator(tchunk_iter);
	return -1;
}

/* Return:
 *     1 = if the chunk before the move was not the last chunk and the move to
 *         the next chunk was successful;
 *     0 = if the chunk before the move was the last chunk and so the move to
 *         the next chunk was not possible;
 *   < 0 = if error
 * Typical use:
 *     while (ret = _next_tchunk(tchunk_iter, tchunk_vps)) {
 *         if (ret < 0) {
 *             an error occured
 *         }
 *         handle 'tchunk_vps'
 *     }
 */
int _next_tchunk(TChunkIterator *tchunk_iter)
{
	const AllTChunks *all_tchunks = tchunk_iter->all_tchunks;
	tchunk_iter->tchunk_rank++;
	if (tchunk_iter->tchunk_rank == all_tchunks->total_num_tchunks)
		return 0;
	const H5DSetDescriptor *h5dset = all_tchunks->h5dset;
	tchunk_iter->moved_along = tchunk_iter->tchunk_rank == 0 ?
					h5dset->ndim :
					next_midx(h5dset->ndim,
						  all_tchunks->num_tchunks,
						  tchunk_iter->tchunk_midx_buf);
	update_TChunkViewports(h5dset,
			       tchunk_iter->tchunk_midx_buf,
			       tchunk_iter->moved_along,
			       all_tchunks->index,
			       all_tchunks->breakpoint_bufs,
			       all_tchunks->tchunkidx_bufs,
			       &tchunk_iter->tchunk_vps);
	return 1;
}

void _print_tchunk_info(const TChunkIterator *tchunk_iter)
{
	const AllTChunks *all_tchunks = tchunk_iter->all_tchunks;
	Rprintf("processing chunk %lld/%lld: [",
		tchunk_iter->tchunk_rank + 1, all_tchunks->total_num_tchunks);
	int ndim = all_tchunks->h5dset->ndim;
	for (int along = 0; along < ndim; along++) {
		size_t i = tchunk_iter->tchunk_midx_buf[along] + 1;
		if (along != 0)
			Rprintf(", ");
		Rprintf("%lu/%lu", i, all_tchunks->num_tchunks[along]);
	}
	Rprintf("] -- <<");
	int along, h5along;
	for (along = 0, h5along = ndim - 1; along < ndim; along++, h5along--) {
		size_t i = tchunk_iter->tchunk_midx_buf[along];
		SEXP start = GET_LIST_ELT(all_tchunks->index, along);
		long long int tchunkidx;
		if (start != R_NilValue) {
			const LLongAEAE *tchunkidx_bufs =
						all_tchunks->tchunkidx_bufs;
			tchunkidx = tchunkidx_bufs->elts[along]->elts[i];
		} else {
			tchunkidx = (long long int) i;
		}
		if (along != 0)
			Rprintf(", ");
		Rprintf("#%lld=%llu:%llu", tchunkidx + 1,
			tchunk_iter->tchunk_vps.h5chunk_vp.h5off[h5along] + 1,
			tchunk_iter->tchunk_vps.h5chunk_vp.h5off[h5along] +
			    tchunk_iter->tchunk_vps.h5chunk_vp.h5dim[h5along]);
	}
	Rprintf(">>\n");
	return;
}


/****************************************************************************
 * ChunkDataBuffer utilities
 */

void _destroy_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf)
{
	if (chunk_data_buf->data_space_id != -1)
		H5Sclose(chunk_data_buf->data_space_id);
	if (chunk_data_buf->data != NULL)
		free(chunk_data_buf->data);
	if (chunk_data_buf->data_vp.h5off != NULL)
		free(chunk_data_buf->data_vp.h5off);
	if (chunk_data_buf->compressed_data != NULL)
		free(chunk_data_buf->compressed_data);
	return;
}

int _init_ChunkDataBuffer(ChunkDataBuffer *chunk_data_buf,
		const H5DSetDescriptor *h5dset, int use_Rtype)
{
	size_t data_length, data_type_size;
	hid_t data_type_id;
	int h5along;
	const H5TypeDescriptor *h5type;

	if (h5dset->h5chunkdim == NULL) {
		PRINT_TO_ERRMSG_BUF("'h5dset->h5chunkdim' is NULL");
		return -1;
	}

	/* Initialize ChunkDataBuffer struct members that control
	   what _destroy_ChunkDataBuffer() needs to free or close. */
	chunk_data_buf->data_space_id = -1;
	chunk_data_buf->data = NULL;
	chunk_data_buf->data_vp.h5off = NULL;
	chunk_data_buf->compressed_data = NULL;

	/* Set struct member 'data_length'. */
	data_length = 1;
	for (h5along = 0; h5along < h5dset->ndim; h5along++)
		data_length *= h5dset->h5chunkdim[h5along];
	chunk_data_buf->data_length = data_length;

	/* Set struct members 'data_type_id' and 'data_type_size'. */
	h5type = h5dset->h5type;
	if (h5type->h5class == H5T_STRING) {
		data_type_id = h5type->h5type_id;
		data_type_size = h5type->h5type_size;
	} else if (use_Rtype ||
		   h5type->native_type_size >= h5type->Rtype_size)
	{
		/* Copying data from 'chunk_data_buf' to final R array will
		   require NO type casting. */
		data_type_id = h5type->native_type_id_for_Rtype;
		data_type_size = h5type->Rtype_size;
	} else {
		/* Copying data from 'chunk_data_buf' to final R array will
		   require type casting. */
		data_type_id = h5type->native_type_id;
		data_type_size = h5type->native_type_size;
	}
	chunk_data_buf->data_type_id = data_type_id;
	chunk_data_buf->data_type_size = data_type_size;

	/* Set struct member 'data_size'. */
	chunk_data_buf->data_size = data_length * data_type_size;
	return 0;
}

int _load_chunk(const H5DSetDescriptor *h5dset,
		const TChunkViewports *tchunk_vps,
		ChunkDataBuffer *chunk_data_buf,
		int use_H5Dread_chunk)
{
	hid_t data_space_id;
	int ret;

	if (chunk_data_buf->data == NULL) {
		chunk_data_buf->data = malloc(chunk_data_buf->data_size);
		if (chunk_data_buf->data == NULL) {
			PRINT_TO_ERRMSG_BUF("failed to allocate memory "
					    "for 'chunk_data_buf->data'");
			return -1;
		}
	}
	if (!use_H5Dread_chunk) {
		if (chunk_data_buf->data_space_id == -1) {
			data_space_id = H5Screate_simple(h5dset->ndim,
							 h5dset->h5chunkdim,
							 NULL);
			if (data_space_id < 0) {
				PRINT_TO_ERRMSG_BUF("H5Screate_simple() "
						    "returned an error");
				return -1;
			}
			chunk_data_buf->data_space_id = data_space_id;
		}
		if (chunk_data_buf->data_vp.h5off == NULL) {
			chunk_data_buf->data_vp.h5off =
				_alloc_hsize_t_buf(h5dset->ndim, 1,
					"'chunk_data_buf->data_vp.h5off'");
			if (chunk_data_buf->data_vp.h5off == NULL)
				return -1;
		}
		chunk_data_buf->data_vp.h5dim = tchunk_vps->h5chunk_vp.h5dim;
		ret = _read_H5Viewport(h5dset,
				&tchunk_vps->h5chunk_vp,
				chunk_data_buf->data_type_id,
				chunk_data_buf->data_space_id,
				chunk_data_buf->data,
				&chunk_data_buf->data_vp);
	} else {
		/* Experimental! */
		if (chunk_data_buf->compressed_data == NULL) {
			chunk_data_buf->compressed_data =
				malloc(chunk_data_buf->data_size +
				       CHUNK_COMPRESSION_OVERHEAD);
			if (chunk_data_buf->compressed_data == NULL) {
				PRINT_TO_ERRMSG_BUF(
					"failed to allocate memory for "
					"'chunk_data_buf->compressed_data'");
				return -1;
			}
		}
		ret = read_h5chunk(h5dset->dset_id,
				   &tchunk_vps->h5chunk_vp,
				   chunk_data_buf);
	}
	return ret;
}

int _reclaim_vlen_bufs(ChunkDataBuffer *chunk_data_buf)
{
	int ret;

	ret = H5Dvlen_reclaim(chunk_data_buf->data_type_id,
			      chunk_data_buf->data_space_id,
			      H5P_DEFAULT, chunk_data_buf->data);
	if (ret < 0)
		PRINT_TO_ERRMSG_BUF("H5Dvlen_reclaim() returned an error");
	return ret;
}

