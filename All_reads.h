#ifndef ALL_READS
#define ALL_READS

#include <stdint.h>

typedef struct {
	uint64_t qns;
	uint32_t qe, tn, ts, te;
	// uint32_t ml:31, rev:1;
	uint32_t cc:30, ml:1, rev:1;
	uint32_t bl:31, del:1;
	uint8_t el;
	uint8_t no_l_indel;
} ma_hit_t;

typedef struct {
	ma_hit_t* buffer;
    uint32_t size;
    uint32_t length;
	uint8_t is_fully_corrected;
	uint8_t is_abnormal;
} ma_hit_t_alloc;

typedef struct
{
    /**[0-1] bits are type:**/
    /**[2-31] bits are length**/
    uint32_t* record;
    uint32_t length;
	uint32_t size;

    char* lost_base;
    uint32_t lost_base_length;
	uint32_t lost_base_size;
	uint32_t new_length;
} Compressed_Cigar_record;


typedef struct
{
	uint64_t** N_site;
	///uint8_t* read;
	char* name;

	uint8_t** read_sperate;
	uint64_t* read_length;
	uint64_t* read_size;
	uint8_t* trio_flag;
    /*
    KJ: 
    in each dirty_reads element: 
        top 2 bits stores the latest round the read was corrected in
        rest of the bits store whether the read was corrected in i-th round
        [10]00 0101 --> latest corrected in round 2, was dirty in rounds 0 and 2 
    */
	uint8_t* dirty_reads;
    uint8_t** rsc;
    uint8_t round;//KJ: TODO: only works up to 8 rounds

	///seq start pos in uint8_t* read
	///do not need it
	///uint64_t* index;
	uint64_t index_size;

    ///name start pos in char* name
	uint64_t* name_index;
	uint64_t name_index_size;
	uint64_t total_reads;
	uint64_t total_reads0;
	uint64_t total_reads_bases;
	uint64_t total_name_length;

	Compressed_Cigar_record* cigars; 
	Compressed_Cigar_record* second_round_cigar;

    ma_hit_t_alloc* paf;
    ma_hit_t_alloc* reverse_paf;

    ///kvec_t_u64_warp* pb_regions;
} All_reads;

#endif //ALL_READS
