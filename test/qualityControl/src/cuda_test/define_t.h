/*start*/
#define TOTAL_SIZE 	70000
#define BLOCKS_NUM	4	//should be divisble by TOTAL_SIZE
#define BLOCK_SIZE	(TOTAL_SIZE/BLOCKS_NUM)

#define Y_ROWS 12
#define Y_COLS TOTAL_SIZE

#define RANDOM_MAX 10
#define RANDOM_MIN 5

#define BLOCK_WIDTH 256
#define TILE_WIDTH 16

#define CUSTOM_VAL 512 //should be a multiple of tile_wdth

typedef long int int_t;

#define CPU_COMPUTE

/*end*/
