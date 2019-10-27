#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <stddef.h>

#ifndef NATIVE
#pragma pack(1)

static void bla()
{
    fputs("called an uninitialized function pointer\n", stderr);
    exit(EXIT_FAILURE);
}

#define SYMBOLT(x, T) T (*p##x)() = (void*)bla;
#include "symbols.inc"
#undef SYMBOLT
#endif

/* .h4m (HVQM4 1.3/1.5) audio decoder 0.4 by flacs/hcs */

//#define VERBOSE_PRINT

/* big endian */

static uint8_t get8(FILE *infile)
{
    uint8_t buf[1];
    if (1 != fread(buf, 1, 1, infile))
    {
        fprintf(stderr, "read error at 0x%lx\n", (unsigned long)ftell(infile));
        exit(EXIT_FAILURE);
    }

    return buf[0];
}

static uint16_t read16(void const * buf)
{
    uint32_t v = 0;
    for (int i = 0; i < 2; i++)
    {
        v <<= 8;
        v |= ((uint8_t const *)buf)[i];
    }
    return v;
}

static uint16_t get16(FILE * infile)
{
    uint8_t buf[2];
    if (2 != fread(buf, 1, 2, infile))
    {
        fprintf(stderr, "read error at 0x%lx\n", (unsigned long)ftell(infile));
        exit(EXIT_FAILURE);
    }
    return read16(buf);
}

static uint32_t read32(void const * buf)
{
    uint32_t v = 0;
    for (int i = 0; i < 4; i++)
    {
        v <<= 8;
        v |= ((uint8_t const*)buf)[i];
    }
    return v;
}

static uint32_t get32(FILE * infile)
{
    uint8_t buf[4];
    if (4 != fread(buf, 1, 4, infile))
    {
        fprintf(stderr, "read error at 0x%lx\n", (unsigned long)ftell(infile));
        exit(EXIT_FAILURE);
    }
    return read32(buf);
}

static void expect32(uint32_t expected, FILE * infile)
{
    uint32_t v = get32(infile);
    if (v != expected)
    {
        fprintf(stderr, "expected 0x%08"PRIx32" at 0x%lx, got 0x%08"PRIx32"\n",
            expected, ftell(infile)-4, v);
        exit(EXIT_FAILURE);
    }
}

static void expect32_imm(uint32_t expected, uint32_t actual, unsigned long offset)
{
    if (expected != actual)
    {
        fprintf(stderr, "expected 0x%08lx to be 0x%08"PRIx32", got 0x%08"PRIx32"\n",
            offset, expected, actual);
        exit(EXIT_FAILURE);
    }
}

static void expect32_text(uint32_t expected, uint32_t actual, char const *name)
{
    if (expected != actual)
    {
        fprintf(stderr, "expected %s to be 0x%08"PRIx32", got 0x%08"PRIx32"\n",
            name, expected, actual);
        exit(EXIT_FAILURE);
    }
}

static void seek_past(uint32_t offset, FILE * infile)
{
    if (-1 == fseek(infile, offset, SEEK_CUR))
    {
        fprintf(stderr, "seek by 0x%x failed\n", offset);
        exit(EXIT_FAILURE);
    }
}

/* audio decode */

struct audio_state
{
    struct
    {
        int16_t hist;
        int8_t idx;
    } *ch;
};

const int32_t IMA_Steps[89] =
{
    7, 8, 9, 10, 11, 12, 13, 14,
    16, 17, 19, 21, 23, 25, 28, 31,
    34, 37, 41, 45, 50, 55, 60, 66,
    73, 80, 88, 97, 107, 118, 130, 143,
    157, 173, 190, 209, 230, 253, 279, 307,
    337, 371, 408, 449, 494, 544, 598, 658,
    724, 796, 876, 963, 1060, 1166, 1282, 1411,
    1552, 1707, 1878, 2066, 2272, 2499, 2749, 3024,
    3327, 3660, 4026, 4428, 4871, 5358, 5894, 6484,
    7132, 7845, 8630, 9493, 10442, 11487, 12635, 13899,
    15289, 16818, 18500, 20350, 22385, 24623, 27086, 29794,
    32767

};

const int8_t IMA_IndexTable[16] = 

{
    -1, -1, -1, -1, 2, 4, 6, 8,
    -1, -1, -1, -1, 2, 4, 6, 8 
};

static int16_t clamp16(int32_t v)
{
    if (v > INT16_MAX) return INT16_MAX;
    else if (v < INT16_MIN) return INT16_MIN;
    return v;
}

static void decode_audio(struct audio_state *state, int first_aud, uint32_t sample_count, FILE *infile, FILE *outfile, int channels)
{
    int16_t * samples = samples = malloc(sample_count*sizeof(int16_t)*channels);
    uint32_t t = 0, i = 0;

    if (first_aud)
    {
        for (int c = channels - 1; c >= 0; c--)
        {
            state->ch[c].hist = get8(infile);
            state->ch[c].hist <<= 8;
            uint8_t b = get8(infile);
            state->ch[c].hist |= (b & 0x80);
            state->ch[c].idx = b & 0x7f;

            if (state->ch[c].idx > 88)
            {
                fprintf(stderr, "invalid step index (%d) at 0x%lx\n", state->ch[c].idx, ftell(infile));
                exit(EXIT_FAILURE);
            }
        }

        for (int c = 0; c < channels; c++)
        {
            samples[t++] = state->ch[c].hist;
        }
        i ++;
    }

    uint8_t b;
    int bitsleft = 0;
    for (; i < sample_count; i++)
    {
        for (int c = channels-1; c >= 0; c--)
        {
            if (bitsleft == 0)
            {
                b = get8(infile);
                bitsleft = 8;
            }

            int32_t ima_step = IMA_Steps[state->ch[c].idx];
            int32_t ima_delta = ima_step >> 3;
            if (b & 0x10) ima_delta += ima_step >> 2;
            if (b & 0x20) ima_delta += ima_step >> 1;
            if (b & 0x40) ima_delta += ima_step;

            if (b & 0x80) state->ch[c].hist = clamp16(state->ch[c].hist - ima_delta);
            else state->ch[c].hist = clamp16(state->ch[c].hist + ima_delta);
            state->ch[c].idx += IMA_IndexTable[(b&0xf0)>>4];
            if (state->ch[c].idx > 88) state->ch[c].idx = 88;
            if (state->ch[c].idx < 0) state->ch[c].idx = 0;

            b <<= 4;
            bitsleft -= 4;
        }

        for (int c = 0; c < channels; c++)
        {
            samples[t++] = state->ch[c].hist;
        }
    }

    if (channels)
    {
        if (sample_count != fwrite(samples, sizeof(int16_t)*channels, sample_count, outfile))
        {
            perror("fwrite");
            exit(EXIT_FAILURE);
        }
    }

    free(samples);
}

/* video stuff */

static int32_t divTable[0x10];
static int32_t mcdivTable[0x200];

static void init_global_constants(void)
{
    divTable[0] = 0;
    mcdivTable[0] = 0;
    for (int i = 1; i < 0x10; ++i)
        divTable[i] = 0x1000 / (i * 16) * 16;
    for (int i = 1; i < 0x200; ++i)
        mcdivTable[i] = 0x1000 / i;
}

static void HVQM4InitDecoder(void)
{
    init_global_constants();
}

// 4x4 block of single value
static void dcBlock(uint8_t *dst, uint32_t stride, uint8_t value)
{
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            dst[y * stride + x] = value;
}

static uint8_t saturate(int32_t x)
{
    return x < 0 ? 0 : x > 0xFF ? 0xFF : x;
}

static uint8_t sat_mean8(uint32_t u)
{
    return saturate((u + 4) / 8);
}

static void WeightImBlock(uint8_t *dst, uint32_t stride, uint8_t value, uint8_t top, uint8_t bottom, uint8_t left, uint8_t right)
{
    /*

    +---+---+---+
    |   | T |   |
    +---+---+---+
    | L | D | R |
    +---+---+---+
    |   | B |   |
    +---+---+---+

     */
    int32_t tmb = top - bottom;
    int32_t lmr = left - right;
    int32_t vph = tmb + lmr;
    int32_t vmh = tmb - lmr;

    int32_t v2 = value * 2;
    int32_t v8 = value * 8;

    int32_t tpl = (top    + left ) - v2;
    int32_t tpr = (top    + right) - v2;
    int32_t bpr = (bottom + right) - v2;
    int32_t bpl = (bottom + left ) - v2;

    int32_t tml = top    - left;
    int32_t tmr = top    - right;
    int32_t bmr = bottom - right;
    int32_t bml = bottom - left;

    // V:
    // 6  8  8 6
    // 8 10 10 8
    // 8 10 10 8
    // 6  8  8 6
    //
    // T:
    //  2  2  2  2
    //  0  0  0  0
    // -1 -1 -1 -1
    // -1 -1 -1 -1
    //
    // B/L/R: like T but rotated accordingly

    // (6*V + 2*T - B + 2*L -   R + 4) / 8
    // (8*V + 2*T - B       -   R + 4) / 8
    // (8*V + 2*T - B -   L       + 4) / 8
    // (6*V + 2*T - B -   L + 2*R + 4) / 8

    dst[0] = sat_mean8(vph + tpl + v8);
    dst[1] = sat_mean8(vph + tml + v8);
    dst[2] = sat_mean8(vmh + tmr + v8);
    dst[3] = sat_mean8(vmh + tpr + v8);

    dst += stride;

    // ( 8*V - B + 2*L -   R + 4) / 8
    // (10*V - B       -   R + 4) / 8
    // (10*V - B -   L       + 4) / 8
    // ( 8*V - B -   L + 2*R + 4) / 8

    dst[0] = sat_mean8(v8 + vph - tml);
    dst[1] = sat_mean8(v8 - bpr      );
    dst[2] = sat_mean8(v8 - bpl      );
    dst[3] = sat_mean8(v8 + vmh - tmr);

    dst += stride;

    // ( 8*V - T + 2*L - R + 4) / 8
    // (10*V - T       - R + 4) / 8
    // (10*V - T - L

    dst[0] = sat_mean8(v8 - vmh - bml);
    dst[1] = sat_mean8(v8 - tpr      );
    dst[2] = sat_mean8(v8 - tpl      );
    dst[3] = sat_mean8(v8 - vph - bmr);

    dst += stride;

    dst[0] = sat_mean8(v8 - vmh + bpl);
    dst[1] = sat_mean8(v8 - vmh + bml);
    dst[2] = sat_mean8(v8 - vph + bmr);
    dst[3] = sat_mean8(v8 - vph + bpr);
}

typedef struct
{
    uint32_t pos;
    int32_t root;
    uint32_t array[2][0x200];
} Tree;

typedef struct
{
    void const *ptr;   // 0-3
    uint32_t size; // 4-7
    uint32_t value;    // 8-B
    int32_t bit;       // C-F
} BitBuffer;
#ifndef NATIVE
_Static_assert(sizeof(BitBuffer) == 0x10, "sizeof(BitBuffer) is incorrect");
#endif

typedef struct
{
    BitBuffer buf; // 0-F
    Tree *tree;    // 0x10-0x13
} BitBufferWithTree;
#ifndef NATIVE
_Static_assert(sizeof(BitBufferWithTree) == 0x14, "sizeof(BitBufferWithTree) is incorrect");
#endif

typedef struct
{
    uint8_t value;
    uint8_t type;
} BlockData;

typedef struct
{
    // size: 0x38 (56)
    BlockData *border; // 0-3 beginning of the plane including the border
    BlockData *payload; // 4-7 beginning of the non-border plane data
    uint16_t h_blocks; // 8-9
    uint16_t v_blocks; // A-B
    uint16_t h_blocks_safe; // C-D
    uint16_t v_blocks_safe; // E-F
    uint16_t some_half_array[4]; // 10-17
    uint32_t some_word_array[4]; // 18-27
    uint16_t width_in_samples; // 28-29
    uint16_t height_in_samples; // 2A-2B
    uint32_t size_in_samples; // 2C-2F
    uint8_t width_shift; // 30
    uint8_t height_shift; // 31
    uint8_t h_samp_per_block; // 32
    uint8_t v_samp_per_block; // 33
    uint8_t block_size_in_samples; // 34
    uint8_t padding[3]; // 35-37
} HVQPlaneDesc;
#ifndef NATIVE
_Static_assert(sizeof(HVQPlaneDesc) == 0x38, "sizeof(HVQPlaneDesc) is incorrect");
#endif

typedef struct
{
    // size: 0x6CD8 (plane data comes right after)
    HVQPlaneDesc planes[3]; // 0x00 - 0xA8
    Tree trees[6];
    BitBufferWithTree symbBuff[3];
    BitBufferWithTree dc_rle[3];
    BitBufferWithTree bufTree0[3];
    BitBufferWithTree bufTree1[2];
    BitBufferWithTree bufTree2[2];
    BitBuffer buf0[3];
    BitBufferWithTree bufTree3[2];
    BitBufferWithTree bufTree4[2]; // 0x6234-0x625B
    uint16_t h_nest_size;
    uint16_t v_nest_size;
    uint8_t is_landscape; // FIXME: check what happens for square video
    uint8_t nest_data[70 * 38]; // 0x6261
    uint8_t padding[3];
    uint32_t boundB; // 0x6CC8
    uint32_t boundA; // 0x6CCC
    uint8_t unk_shift; // 0x6CD0
    uint8_t unk6CD1; // 0x6CD1
    uint8_t unk6CD2[2]; // 0x6CD2
    uint8_t unk6CD4[2]; // 0x6CD4
    uint8_t maybe_padding[2]; // 0x6CD6-0x6CD7
} VideoState;
#ifndef NATIVE
_Static_assert(sizeof(VideoState) == 0x6CD8, "sizeof(VideoState) is incorrect");
#endif

typedef struct
{
    VideoState *state;
    uint16_t width;
    uint16_t height;
    uint8_t h_samp;
    uint8_t v_samp;
    uint16_t unkA; // unused?
} SeqObj;

typedef struct
{
    SeqObj seqobj;
    void *past;
    void *present;
    void *future;
} Player;

typedef struct
{
    uint16_t hres;
    uint16_t vres;
    uint8_t h_samp;
    uint8_t v_samp;
    uint8_t video_mode;
} VideoInfo;

// copy into 4x4 block
static void OrgBlock(VideoState *state, uint8_t *dst, uint32_t dst_stride, uint32_t plane_idx)
{
    BitBuffer *buf = &state->buf0[plane_idx];
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            dst[y * dst_stride + x] = *(uint8_t*)buf->ptr++;
}

// can use FFmpeg's get_bits1()/get_bits(..., 8) for this
static int16_t getBit(BitBuffer *buf)
{
    int32_t bit = buf->bit;
    if (bit < 0)
    {
        buf->value = read32(buf->ptr);
        buf->ptr += 4;
        bit = 31;
    }
    buf->bit = bit - 1;
    return (buf->value >> bit) & 1;
}

static uint8_t getByte(BitBuffer *buf)
{
    uint32_t value = buf->value;
    int32_t bit = buf->bit;
    if (bit < 7)
    {
        buf->value = read32(buf->ptr);
        buf->ptr += 4;
        value <<= 7 - bit;
        value |= buf->value >> (bit + 25);
        bit += 24;
    }
    else
    {
        value >>= bit - 7;
        bit -= 8;
    }
    buf->bit = bit;
    return value & 0xFF;
}

static uint32_t readTree_signed;
static uint32_t readTree_scale;

static int16_t _readTree(Tree *dst, BitBuffer *src)
{
    if (getBit(src) == 0)
    {
        // leaf node
        uint8_t byte = getByte(src);
        uint32_t symbol = byte;
        if (readTree_signed && byte > 0x7F)
            symbol = (int8_t)byte;
        symbol <<= readTree_scale;
        dst->array[0][byte] = symbol;
        return byte;
    }
    else
    {
        // recurse
        uint32_t pos = dst->pos++;
        // read the 0 side of the tree
        dst->array[0][pos] = (uint32_t)_readTree(dst, src);
        // read the 1 side of the tree
        dst->array[1][pos] = (uint32_t)_readTree(dst, src);
        return (int16_t)pos;
    }
}

static void readTree(BitBufferWithTree *buf, uint32_t isSigned, uint32_t scale)
{
    readTree_signed = isSigned;
    readTree_scale = scale;
    Tree *tree = buf->tree;
    tree->pos = 0x100;
    if (buf->buf.size == 0)
        tree->root = 0;
    else
        tree->root = _readTree(tree, &buf->buf);
}

static uint32_t decodeHuff(BitBufferWithTree *buf)
{
    Tree *tree = buf->tree;
    int32_t pos = tree->root;
    while (pos >= 0x100)
        pos = tree->array[getBit(&buf->buf)][pos];
    return tree->array[0][pos];
}

static int32_t decodeSOvfSym(BitBufferWithTree *buf, int32_t a, int32_t b)
{
    int32_t sum = 0;
    int32_t value;
    do
    {
        value = decodeHuff(buf);
        sum += value;
    } while (value <= a || value >= b);
    return sum;
}

static int32_t decodeUOvfSym(BitBufferWithTree *buf, int32_t cmp_value)
{
    int32_t sum = 0;
    int32_t value;
    do
    {
        value = decodeHuff(buf);
        sum += value;
    } while (value >= cmp_value);
    return sum;
}

static uint32_t GetAotBasis(VideoState *state, uint8_t dst[4][4], int32_t *sum, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    BitBuffer *buf = &state->buf0[plane_idx];
    uint16_t bits = read16(buf->ptr);
    buf->ptr += 2;
    uint32_t x_stride, y_stride;
    uint32_t big = bits & 0x3F;
    uint32_t small = (bits >> 6) & 0x1F;
    if (state->is_landscape)
    {
        nest_data += nest_stride * small + big;
        x_stride =           1 << ((bits >> 11) & 1);
        y_stride = nest_stride << ((bits >> 12) & 1);
    }
    else
    {
        nest_data += nest_stride * big + small;
        x_stride =           1 << ((bits >> 12) & 1);
        y_stride = nest_stride << ((bits >> 11) & 1);
    }
    uint8_t min = *nest_data;
    uint8_t max = *nest_data;
    for (int y = 0; y < 4; ++y)
    {
        for (int x = 0; x < 4; ++x)
        {
            uint8_t nest_value = nest_data[y * y_stride + x * x_stride];
            dst[y][x] = nest_value;
            min = nest_value < min ? nest_value : min;
            max = nest_value > max ? nest_value : max;
        }
    }
    *sum += decodeHuff(&state->bufTree0[plane_idx]);
    int32_t inverse = divTable[max - min];
    if (bits & 0x8000)
        inverse = -inverse;
    int32_t foo = (bits >> 13) & 3;
    return (*sum + foo) * inverse;
}

static uint32_t GetMCAotBasis(VideoState *state, uint8_t dst[4][4], int32_t *sum, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    // the only difference to GetAotBasis() seems to be the ">> 4 & 0xF"
    BitBuffer *buf = &state->buf0[plane_idx];
    uint16_t bits = read16(buf->ptr);
    buf->ptr += 2;
    uint32_t step, stride;
    uint32_t big = bits & 0x3F;
    uint32_t small = (bits >> 6) & 0x1F;
    if (state->is_landscape)
    {
        nest_data += nest_stride * small + big;
        step   =           1 << ((bits >> 11) & 1);
        stride = nest_stride << ((bits >> 12) & 1);
    }
    else
    {
        nest_data += nest_stride * big + small;
        step   =           1 << ((bits >> 12) & 1);
        stride = nest_stride << ((bits >> 11) & 1);
    }
    uint8_t min, max;
    min = max = (nest_data[0] >> 4) & 0xF; // !
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            uint8_t nest_value = (nest_data[i * stride + j * step] >> 4) & 0xF; // !
            dst[i][j] = nest_value;
            min = nest_value < min ? nest_value : min;
            max = nest_value > max ? nest_value : max;
        }
    }
    *sum += decodeHuff(&state->bufTree0[plane_idx]);
    int32_t inverse = divTable[max - min];
    if (bits & 0x8000)
        inverse = -inverse;
    int32_t foo = (bits >> 13) & 3;
    return (*sum + foo) * inverse;
}

static int32_t GetAotSum(VideoState *state, int32_t result[4][4], uint8_t count, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            result[y][x] = 0;
    uint8_t byte_result[4][4];
    int32_t temp = 0;
    for (int k = 0; k < count; ++k)
    {
        uint32_t factor = GetAotBasis(state, byte_result, &temp, nest_data, nest_stride, plane_idx);
        for (int y = 0; y < 4; ++y)
            for (int x = 0; x < 4; ++x)
                result[y][x] += factor * byte_result[y][x];
    }
    int32_t sum = 0;
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            sum += result[y][x];
    return sum >> 4;
}

static int32_t GetMCAotSum(VideoState *state, int32_t result[4][4], uint8_t count, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            result[i][j] = 0;
    uint8_t byte_result[4][4];
    int32_t temp = 0;
    for (int k = 0; k < count; ++k)
    {
        uint32_t factor = GetMCAotBasis(state, byte_result, &temp, nest_data, nest_stride, plane_idx);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                result[i][j] += factor * byte_result[i][j];
    }
    int32_t sum = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            sum += result[i][j];
    return sum >> 4;
}

static void HVQM4InitSeqObj(SeqObj *seqobj, VideoInfo *videoinfo)
{
    seqobj->width = videoinfo->hres;
    seqobj->height = videoinfo->vres;
    // TODO: find better names for these, uv_h_step/uv_v_step?
    seqobj->h_samp = videoinfo->h_samp;
    seqobj->v_samp = videoinfo->v_samp;
}

static uint32_t HVQM4BuffSize(SeqObj *seqobj)
{
    uint32_t h_blocks = seqobj->width / 4;
    uint32_t v_blocks = seqobj->height / 4;
    uint32_t y_blocks = (h_blocks + 2) * (v_blocks + 2);

    uint32_t uv_h_blocks = seqobj->h_samp == 2 ? h_blocks / 2 : h_blocks;
    uint32_t uv_v_blocks = seqobj->v_samp == 2 ? v_blocks / 2 : v_blocks;
    uint32_t uv_blocks = (uv_h_blocks + 2) * (uv_v_blocks + 2);

    uint32_t total = (y_blocks + uv_blocks * 2) * sizeof(uint16_t);
    return sizeof(VideoState) + total;
}

// h_samp/v_samp: pixels per sample
static void setHVQPlaneDesc(SeqObj *seqobj, uint8_t plane_idx, uint8_t h_samp, uint8_t v_samp)
{
    HVQPlaneDesc *plane = &seqobj->state->planes[plane_idx];
    plane->width_shift = h_samp == 2 ? 1 : 0;
    plane->width_in_samples = seqobj->width >> plane->width_shift;
    plane->height_shift = v_samp == 2 ? 1 : 0;
    plane->height_in_samples = seqobj->height >> plane->height_shift;
    plane->size_in_samples = plane->width_in_samples * plane->height_in_samples;
    // pixels per 2x2 block
    plane->h_samp_per_block = 2 >> plane->width_shift;
    plane->v_samp_per_block = 2 >> plane->height_shift;
    plane->block_size_in_samples = plane->h_samp_per_block * plane->v_samp_per_block;
    // number of 4x4 blocks
    plane->h_blocks = seqobj->width / (h_samp * 4);
    plane->v_blocks = seqobj->height / (v_samp * 4);
    // number of 4x4 blocks + border
    plane->h_blocks_safe = plane->h_blocks + 2;
    plane->v_blocks_safe = plane->v_blocks + 2;
    plane->some_half_array[0] = 0;
    plane->some_half_array[1] = plane->h_blocks_safe;
    plane->some_half_array[2] = plane->h_blocks_safe + 1;
    plane->some_half_array[3] = 1;
    plane->some_word_array[0] = 0;
    plane->some_word_array[1] = plane->width_in_samples << 2;
    plane->some_word_array[2] = (plane->width_in_samples << 2) + 4;
    plane->some_word_array[3] = 4;
}

// HACK: assumes 4:2:0
__attribute__((unused))
static void dumpYUV(Player *player, char const *path)
{
    FILE *f = fopen(path, "wb+");
    uint32_t w = player->seqobj.width, h = player->seqobj.height;
    fprintf(f, "P5\n%u %u\n255\n", w, h * 2);
    uint8_t const *p = player->present;
    for (int plane = 0; plane < 2; ++plane)
    {
        for (uint32_t i = 0; i < h; ++i)
        {
            for (uint32_t j = 0; j < w; ++j)
            {
                if (plane == 0 || j < w / 2)
                    fputc(*p++, f);
                else
                    fputc(0, f);
            }
        }
    }
    fclose(f);
}

// HACK: assumes 4:2:0, assumes JPEG color space
static uint8_t clamp255(float f)
{
    return f < 0 ? 0 : f > 255 ? 255 : (uint8_t)f;
}
static void dumpRGB(Player *player, char const *path)
{
    FILE *f = fopen(path, "wb+");
    uint32_t w = player->seqobj.width, h = player->seqobj.height;
    fprintf(f, "P6\n%u %u\n255\n", w, h);
    uint8_t const *yp = player->present;
    uint8_t const *up = yp + w*h;
    uint8_t const *vp = up + w*h/4;
    uint8_t *rgb = malloc(w * h * 3);
    uint8_t *ptr = rgb;
    for (uint32_t i = 0; i < h; ++i)
    {
        for (uint32_t j = 0; j < w; ++j)
        {
            float y = yp[i   * w   + j];
            float u = up[i/2 * w/2 + j/2];
            float v = vp[i/2 * w/2 + j/2];
            *ptr++ = clamp255(y + 1.402f*(v - 128.f));
            *ptr++ = clamp255(y - 0.34414f*(u - 128.f) - 0.71414f*(v - 128.f));
            *ptr++ = clamp255(y + 1.772f*(u - 128.f));
        }
    }
    fwrite(rgb, w*h*3, 1, f);
    fclose(f);
    free(rgb);
}

__attribute__((unused))
static void dumpPlanes(VideoState *state, char const *prefix)
{
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        char path[128];
        snprintf(path, 128, "%s_%c.ppm", prefix, "yuv"[plane_idx]);
        FILE *f = fopen(path, "wb+");
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        fprintf(f, "P5\n%u %u\n255\n", plane->h_blocks_safe, plane->v_blocks_safe);
        uint8_t const *p = (uint8_t const*)plane->border;
        for (int i = 0; i < plane->v_blocks_safe; ++i)
        {
            for (int j = 0; j < plane->h_blocks_safe; ++j)
            {
                fputc(*p++, f);
                ++p;
            }
        }
        fclose(f);
    }
}

static void set_border(BlockData *dst)
{
    dst->value = 0x7F;
    dst->type = 0xFF;
}

static void HVQM4SetBuffer(SeqObj *seqobj, void *workbuff)
{
    VideoState *state = workbuff;
    seqobj->state = state;
    setHVQPlaneDesc(seqobj, 0, 1, 1);
    setHVQPlaneDesc(seqobj, 1, seqobj->h_samp, seqobj->v_samp);
    setHVQPlaneDesc(seqobj, 2, seqobj->h_samp, seqobj->v_samp);

    state->is_landscape = seqobj->width >= seqobj->height;
    if (state->is_landscape)
    {
        state->h_nest_size = 70;
        state->v_nest_size = 38;
    }
    else
    {
        state->h_nest_size = 38;
        state->v_nest_size = 70;
    }

    state->bufTree1[0].tree = &state->trees[3];
    state->bufTree1[1].tree = &state->trees[3];

    state->bufTree2[0].tree = &state->trees[1];
    state->bufTree2[1].tree = &state->trees[1];

    state->symbBuff[0].tree = &state->trees[0];
    state->symbBuff[1].tree = &state->trees[0];
    state->symbBuff[2].tree = &state->trees[0];

    state->dc_rle[0].tree = &state->trees[1]; // reuse!
    state->dc_rle[1].tree = &state->trees[1]; //
    state->dc_rle[2].tree = &state->trees[1]; //

    state->bufTree0[0].tree = &state->trees[2];
    state->bufTree0[1].tree = &state->trees[2];
    state->bufTree0[2].tree = &state->trees[2];

    state->bufTree3[0].tree = &state->trees[4];
    state->bufTree3[1].tree = &state->trees[4];

    state->bufTree4[0].tree = &state->trees[5];
    state->bufTree4[1].tree = &state->trees[5];

    BlockData *plane_data = workbuff + sizeof(VideoState);
    for (int i = 0; i < 3; ++i)
    {
        HVQPlaneDesc *plane = &state->planes[i];
        plane->border = plane_data;
        uint32_t stride = plane->h_blocks_safe;
        // skip top border (stride) and left border (1)
        plane->payload = plane_data + stride + 1;
        plane_data += plane->h_blocks_safe * plane->v_blocks_safe;

        // set horizontal borders
        BlockData *ptr = plane->border;
        for (uint32_t i = plane->h_blocks_safe; i; --i)
        {
            set_border(ptr);
            ++ptr;
        }

        ptr = plane_data;
        for (uint32_t i = plane->h_blocks_safe; i; --i)
        {
            --ptr;
            set_border(ptr);
        }

        // set vertical borders
        ptr = plane->border + stride;
        for (uint32_t i = plane->v_blocks_safe - 2; i; --i)
        {
            set_border(ptr);
            ptr += stride;
        }

        ptr = plane->border + stride * 2 - 1;
        for (uint32_t i = plane->v_blocks_safe - 2; i; --i)
        {
            set_border(ptr);
            ptr += stride;
        }
    }
}

static uint32_t getDeltaDC(VideoState *state, uint32_t plane_idx, uint32_t *rle_state)
{
    if (*rle_state == 0)
    {
        uint32_t delta = decodeSOvfSym(&state->symbBuff[plane_idx], state->boundA, state->boundB);
        // successive zeroes are run-length encoded
        if (delta == 0)
            *rle_state = decodeHuff(&state->dc_rle[plane_idx]);
        return delta;
    }
    else
    {
        --*rle_state;
        return 0;
    }
}

// initialize a bit buffer
static void setCode(BitBuffer *dst, void const *src)
{
    dst->size = read32(src);
    dst->ptr = dst->size ? src + 4 : NULL;
    dst->bit = -1;
}

static void Ipic_BasisNumDec(VideoState *state)
{
    BlockData *y_dst = state->planes[0].payload;
    uint32_t y_v_blocks = state->planes[0].v_blocks;
    BitBufferWithTree *y_tree1 = &state->bufTree1[0];
    BitBufferWithTree *y_tree2 = &state->bufTree2[0];
    uint32_t rle_state = 0;
    while (y_v_blocks--)
    {
        uint32_t y_h_blocks = state->planes[0].h_blocks;
        while (y_h_blocks--)
        {
            if (rle_state)
            {
                y_dst->type = 0;
                ++y_dst;
                --rle_state;
            }
            else
            {
                int16_t x = decodeHuff(y_tree1) & 0xFFFF;
                if (x == 0)
                    rle_state = decodeHuff(y_tree2);
                y_dst->type = x & 0xFF;
                ++y_dst;
            }
        }
        // skip borders
        y_dst += 2;
    }

    BlockData *u_dst = state->planes[1].payload;
    BlockData *v_dst = state->planes[2].payload;
    uint32_t uv_v_blocks = state->planes[1].v_blocks;
    BitBufferWithTree *uv_tree1 = &state->bufTree1[1];
    BitBufferWithTree *uv_tree2 = &state->bufTree2[1];
    rle_state = 0;
    while (uv_v_blocks--)
    {
        uint32_t uv_h_blocks = state->planes[1].h_blocks;
        while (uv_h_blocks--)
        {
            if (rle_state)
            {
                u_dst->type = 0;
                v_dst->type = 0;
                ++u_dst;
                ++v_dst;
                --rle_state;
            }
            else
            {
                int16_t value = decodeHuff(uv_tree1) & 0xFFFF;
                if (value == 0)
                    rle_state = decodeHuff(uv_tree2);
                u_dst->type = (value >> 0) & 0xF;
                v_dst->type = (value >> 4) & 0xF;
                ++u_dst;
                ++v_dst;
            }
        }
        u_dst += 2;
        v_dst += 2;
    }
}

static void IpicDcvDec(VideoState *state)
{
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        uint32_t rle_state = 0;
        uint32_t v_blocks = plane->v_blocks;
        BlockData *curr = plane->payload;
        while (v_blocks--)
        {
            // pointer to previous line
            BlockData const *prev = curr - plane->h_blocks_safe;
            // first prediction on a line is only the previous line's value
            uint8_t value = prev->value;
            for (uint32_t i = 0; i < plane->h_blocks; ++i)
            {
                value += getDeltaDC(state, plane_idx, &rle_state);
                curr->value = value;
                ++curr;
                ++prev;
                // next prediction on this line is the mean of left and top values
                // +---+---+
                // |   | T |
                // +---+---+
                // | L | P |
                // +---+---+
                value = (value + prev->value + 1) / 2;
            }
            // skip right border of this line and left border of next line
            curr += 2;
        }
        //dumpPlanes(state, "filled");
    }
}

static void MakeNest(VideoState *state, uint16_t nest_x, uint16_t nest_y)
{
    HVQPlaneDesc *y_plane = &state->planes[0];
    BlockData const *ptr = y_plane->payload + y_plane->h_blocks_safe * nest_y + nest_x;

    int32_t v_empty, h_empty, v_nest_blocks, h_nest_blocks, v_mirror, h_mirror;

    if (y_plane->h_blocks < state->h_nest_size)
    {
        h_nest_blocks = y_plane->h_blocks;
        h_mirror = state->h_nest_size - y_plane->h_blocks;
        if (h_mirror > y_plane->h_blocks)
            h_mirror = y_plane->h_blocks;
        h_empty = state->h_nest_size - (h_nest_blocks + h_mirror);
    }
    else
    {
        h_nest_blocks = state->h_nest_size;
        h_empty = 0;
        h_mirror = 0;
    }

    if (y_plane->v_blocks < state->v_nest_size)
    {
        v_nest_blocks = y_plane->v_blocks;
        v_mirror = state->v_nest_size - y_plane->v_blocks;
        if (v_mirror > y_plane->v_blocks)
            v_mirror = y_plane->v_blocks;
        v_empty = state->v_nest_size - (v_nest_blocks + v_mirror);
    }
    else
    {
        v_nest_blocks = state->v_nest_size;
        v_empty = 0;
        v_mirror = 0;
    }

    uint8_t *nest = state->nest_data;
    for (int i = 0; i < v_nest_blocks; ++i)
    {
        BlockData const *p = ptr;
        for (int j = 0; j < h_nest_blocks; ++j)
        {
            *nest++ = (p->value >> 4) & 0xF;
            ++p;
        }
        for (int j = 0; j < h_mirror; ++j)
        {
            --p;
            *nest++ = (p->value >> 4) & 0xF;
        }
        for (int j = 0; j < h_empty; ++j)
            *nest++ = 0;
        ptr += y_plane->h_blocks_safe;
    }

    uint8_t const *nest2 = nest - state->h_nest_size;
    for (int i = 0; i < v_mirror; ++i)
    {
        for (int j = 0; j < state->h_nest_size; ++j)
            *nest++ = nest2[j];
        nest2 -= state->h_nest_size;
    }

    for (int i = 0; i < v_empty; ++i)
        for (int j = 0; j < state->h_nest_size; ++j)
            *nest++ = 0;
}

static void _MotionComp_00(uint8_t *dst, uint32_t dst_stride, uint8_t const *src, uint32_t src_stride)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            dst[i * dst_stride + j] = src[i * src_stride + j];
}

static void _MotionComp_01(uint8_t *dst, uint32_t dst_stride, uint8_t const *src, uint32_t src_stride)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            dst[i * dst_stride + j] = (src[i * src_stride + j] + src[(i + 1) * src_stride + j] + 1) >> 1;
}

static void _MotionComp_10(uint8_t *dst, uint32_t dst_stride, uint8_t const *src, uint32_t src_stride)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            dst[i * dst_stride + j] = (src[i * src_stride + j] + src[i * src_stride + j + 1] + 1) >> 1;
}

static void _MotionComp_11(uint8_t *dst, uint32_t dst_stride, uint8_t const *src, uint32_t src_stride)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            dst[i * dst_stride + j] = (src[i * src_stride + j] + src[i * src_stride + j + 1] + src[(i + 1) * src_stride + j] + src[(i + 1) * src_stride + j + 1] + 2) >> 2;
}

static void _MotionComp(void *dst, uint32_t dst_stride, void const *src, uint32_t src_stride, uint32_t unk7, uint32_t unk8)
{
    if (unk8 == 0)
        if (unk7 == 0)
            _MotionComp_00(dst, dst_stride, src, src_stride);
        else
            _MotionComp_10(dst, dst_stride, src, src_stride);
    else
        if (unk7 == 0)
            _MotionComp_01(dst, dst_stride, src, src_stride);
        else
            _MotionComp_11(dst, dst_stride, src, src_stride);
}

typedef struct
{
    uint32_t unk0; // init 0
    uint32_t unk4; // init 0x7F
    void *payload_ptr8;
    void *payload_ptrC;
    void *present;
    void *top;
    void *target;
    void *past;
    void *future;
    uint16_t h_unk24;
    uint16_t padding; // ?
    uint32_t v_unk28;
    uint32_t h_samp_per_block;
    uint32_t stride;
} MCPlane;
#ifndef NATIVE
_Static_assert(sizeof(MCPlane) == 0x34, "sizeof(MCPlane) is wrong");
#endif

static void MotionComp(VideoState *state, MCPlane mcplanes[3], int32_t unk5, int32_t unk6)
{
    for (int i = 0; i < 3; ++i)
    {
        MCPlane *mcplane = &mcplanes[i];
        HVQPlaneDesc *plane = &state->planes[i];
        int32_t foo = unk5 >> plane->width_shift;
        int32_t bar = unk6 >> plane->height_shift;
        void *ptr = mcplane->target + (bar >> 1) * plane->width_in_samples + (foo >> 1);
        for (int j = 0; j < plane->block_size_in_samples; ++j)
        {
            _MotionComp(mcplane->top + plane->some_word_array[j],
                        plane->width_in_samples,
                        ptr + plane->some_word_array[j],
                        plane->width_in_samples,
                        foo & 1, bar & 1);
        }
    }
}

static void IntraAotBlock(VideoState *state, uint8_t *dst, uint32_t stride, uint8_t value, uint8_t block_type, uint32_t plane_idx)
{
    if (block_type == 6)
    {
        OrgBlock(state, dst, stride, plane_idx);
        return;
    }
    int32_t r30 = value << state->unk_shift;
    int32_t result[4][4];
    r30 -= GetAotSum(state, result, block_type, state->nest_data, state->h_nest_size, plane_idx);
    for (int y = 0; y < 4; ++y)
    {
        for (int x = 0; x < 4; ++x)
        {
            uint32_t value = ((result[y][x] + r30) >> state->unk_shift);
            dst[y * stride + x] = saturate(value);
        }
    }
}

static void PrediAotBlock(VideoState *state, uint8_t *dst, uint8_t const *src, uint32_t stride, uint8_t block_type,
                          uint8_t *nest_data, uint32_t h_nest_size, uint32_t plane_idx, uint32_t foo, uint32_t bar)
{
    int32_t result[4][4];
    --block_type;
    uint32_t r20 = GetMCAotSum(state, result, block_type, nest_data, h_nest_size, plane_idx);

    uint8_t mdst[4][4];
    uint32_t const dst_stride = 4;
    _MotionComp(mdst, dst_stride, src, stride, foo, bar);
    int32_t mean = 8;
    for (int y = 0; y < 4; ++y)
        for (int x = 0; x < 4; ++x)
            mean += mdst[y][x];
    mean /= 16;
    int32_t diff[4][4];
    int32_t min, max;
    min = max = mdst[0][0] - mean;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int32_t value = diff[i][j] = mdst[i][j] - mean;
            min = value < min ? value : min;
            max = value > max ? value : max;
        }
    }
    uint32_t r28 = (decodeSOvfSym(&state->symbBuff[plane_idx], state->boundA, state->boundB) >> state->unk6CD1 << state->unk_shift) - r20;
    uint32_t bla = (decodeSOvfSym(&state->symbBuff[plane_idx], state->boundA, state->boundB) >> state->unk6CD1);
    uint32_t r18 = bla * mcdivTable[max - min];
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            result[i][j] += r28 + diff[i][j] * r18;

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            uint32_t value = (result[i][j] >> state->unk_shift) + mdst[i][j];
            dst[i * stride + j] = saturate(value);
        }
    }
}

typedef struct
{
    uint32_t plane_idx;
    BlockData const *line_prev;
    BlockData const *line_curr;
    BlockData const *line_next;
    BlockData next;
    BlockData curr;
    uint8_t value_prev;
} StackState;

static void IpicBlockDec(VideoState *state, uint8_t *dst, uint32_t stride, StackState *stack_state)
{
    if (stack_state->curr.type == 0)
    {
        uint8_t top    = stack_state->line_prev->type & 0x77 ? stack_state->curr.value : stack_state->line_prev->value;
        uint8_t bottom = stack_state->line_next->type & 0x77 ? stack_state->curr.value : stack_state->line_next->value;
        uint8_t right  = stack_state->next.type       & 0x77 ? stack_state->curr.value : stack_state->next.value;
        uint8_t left   = stack_state->value_prev;
        WeightImBlock(dst, stride, stack_state->curr.value, top, bottom, left, right);
        stack_state->value_prev = stack_state->curr.value;
    }
    else if (stack_state->curr.type == 8)
    {
        dcBlock(dst, stride, stack_state->curr.value);
        stack_state->value_prev = stack_state->curr.value;
    }
    else
    {
        IntraAotBlock(state, dst, stride, stack_state->curr.value, stack_state->curr.type, stack_state->plane_idx);
        stack_state->value_prev = stack_state->next.value; // ?!
    }
    // next block
    ++stack_state->line_prev;
    ++stack_state->line_next;
}

static void IpicLineDec(VideoState *state, uint8_t *dst, uint32_t stride, StackState *stack_state, uint16_t h_blocks)
{
    stack_state->next = stack_state->line_curr[0];
    stack_state->value_prev = stack_state->line_curr[0].value;

    while (--h_blocks > 0)
    {
        stack_state->curr = stack_state->next;
        ++stack_state->line_curr;
        stack_state->next = stack_state->line_curr[0];
        IpicBlockDec(state, dst, stride, stack_state);
        // next block on same line
        dst += 4;
    }

    stack_state->curr = stack_state->next;
    IpicBlockDec(state, dst, stride, stack_state);

    // skip current, right border on same line, and left border on next line
    stack_state->line_curr += 3;

    // these have already been advanced to the right border in IpicBlockDec
    stack_state->line_prev += 2;
    stack_state->line_next += 2;
}

static void IpicPlaneDec(VideoState *state, int plane_idx, uint8_t *dst)
{
    HVQPlaneDesc *plane = &state->planes[plane_idx];
    StackState stack_state;
    stack_state.plane_idx = plane_idx;
    stack_state.line_prev = plane->payload;
    stack_state.line_curr = plane->payload;
    stack_state.line_next = plane->payload + plane->h_blocks_safe;
    int16_t v_blocks = plane->v_blocks;
    // first line
    if (v_blocks > 0)
    {
        IpicLineDec(state, dst, plane->width_in_samples, &stack_state, plane->h_blocks);
        // blocks are 4x4 so advance dst by 4 lines
        dst += plane->width_in_samples * 4;
        --v_blocks;
    }
    // middle lines
    stack_state.line_prev = plane->payload;
    while (v_blocks > 1)
    {
        IpicLineDec(state, dst, plane->width_in_samples, &stack_state, plane->h_blocks);
        dst += plane->width_in_samples * 4;
        --v_blocks;
    }
    // last line
    if (v_blocks > 0)
    {
        stack_state.line_next = stack_state.line_curr;
        IpicLineDec(state, dst, plane->width_in_samples, &stack_state, plane->h_blocks);
    }
}

static void initMCHandler(VideoState *state, MCPlane mcplanes[3], void *present, void *past, void *future)
{
    for (int i = 0; i < 3; ++i)
    {
        MCPlane *mcplane = &mcplanes[i];
        HVQPlaneDesc *plane = &state->planes[i];
        mcplane->unk0 = 0;
        mcplane->unk4 = 0x7F;
        mcplane->present = present;
        mcplane->past    = past;
        mcplane->future  = future;
        mcplane->payload_ptr8 = plane->payload;
        mcplane->payload_ptrC = plane->payload;
        mcplane->h_unk24 = 8 >> plane->width_shift;
        mcplane->v_unk28 = plane->width_in_samples * (8 >> plane->height_shift);
        mcplane->h_samp_per_block = plane->h_samp_per_block;
        mcplane->stride = plane->h_blocks_safe * plane->v_samp_per_block;

        present += plane->size_in_samples;
        past    += plane->size_in_samples;
        future  += plane->size_in_samples;
    }
}

static void initMCBproc(BitBufferWithTree *buftree, uint32_t *proc)
{
    if (buftree->buf.ptr)
    {
        proc[0] = getBit(&buftree->buf);
        proc[1] = decodeUOvfSym(buftree, 0xFF);
    }
}

static void initMCBtype(BitBufferWithTree *buftree, uint32_t *type)
{
    if (buftree->buf.ptr)
    {
        type[0] = (getBit(&buftree->buf) << 1) | getBit(&buftree->buf);
        type[1] = decodeUOvfSym(buftree, 0xFF);
    }
}

static void setMCTop(MCPlane mcplanes[3])
{
    for (int i = 0; i < 3; ++i)
        mcplanes[i].top = mcplanes[i].present;
}

static uint32_t mcbtypetrans[2][3] = {
    { 1, 2, 0 },
    { 2, 0, 1 },
};

static uint32_t getMCBtype(BitBufferWithTree *buftree, uint32_t *type)
{
    if (type[1] == 0)
    {
        uint32_t bit = getBit(&buftree->buf);
        type[0] = mcbtypetrans[bit][type[0]];
        type[1] = decodeUOvfSym(buftree, 0xFF);
    }
    --type[1];
    return type[0];
}

static uint32_t getMCBproc(BitBufferWithTree *buftree, uint32_t *proc)
{
    if (proc[1] == 0)
    {
        proc[0] ^= 1;
        proc[1] = decodeUOvfSym(buftree, 0xFF);
    }
    --proc[1];
    return proc[0];
}

static void setMCNextBlk(MCPlane mcplanes[3])
{
    for (int i = 0; i < 3; ++i)
    {
        mcplanes[i].top += mcplanes[i].h_unk24;
        mcplanes[i].payload_ptr8 += mcplanes[i].h_samp_per_block * 2;
    }
}

static void setMCDownBlk(MCPlane mcplanes[3])
{
    for (int i = 0; i < 3; ++i)
    {
        MCPlane *mcplane = &mcplanes[i];
        mcplane->present += mcplane->v_unk28;
        void *ptr = mcplane->payload_ptrC + mcplane->stride * 2;
        mcplane->payload_ptr8 = ptr;
        mcplane->payload_ptrC = ptr;
    }
}

static void decode_PB_dc(VideoState *state, MCPlane mcplanes[3])
{
    for (int i = 0; i < 3; ++i)
    {
        HVQPlaneDesc *plane = &state->planes[i];
        MCPlane *mcplane = &mcplanes[i];
        for (int j = 0; j < plane->block_size_in_samples; ++j)
        {
            mcplane->unk4 += decodeSOvfSym(&state->symbBuff[i], state->boundA, state->boundB);
            uint8_t *payload = mcplane->payload_ptr8;
            payload[plane->some_half_array[j] * 2] = mcplane->unk4;
        }
    }
}

static void reset_PB_dc(MCPlane mcplanes[3])
{
    for (int i = 0; i < 3; ++i)
        mcplanes[i].unk4 = 0x7F;
}

static void decode_PB_cc(VideoState *state, MCPlane mcplanes[3], uint32_t unk5, uint32_t unk6)
{
    uint32_t r30 = (unk6 << 5) | (unk5 << 4);
    if (unk5 == 1)
    {
        for (int i = 0; i < 3; ++i)
        {
            uint8_t *payload = mcplanes[i].payload_ptr8;
            HVQPlaneDesc *plane = &state->planes[i];
            for (int j = 0; j < plane->block_size_in_samples; ++j)
                payload[plane->some_half_array[j] * 2 + 1] = r30;
        }
        return;
    }
    else
    {
        HVQPlaneDesc *planeY = &state->planes[0];
        MCPlane *mcplaneY = &mcplanes[0];
        for (int i = 0; i < planeY->block_size_in_samples; ++i)
        {
            uint8_t *ptr = mcplaneY->payload_ptr8;
            if (mcplaneY->unk0)
            {
                ptr[planeY->some_half_array[i] * 2 + 1] = r30;
                --mcplaneY->unk0;
            }
            else
            {
                int16_t huff = decodeHuff(&state->bufTree1[0]);
                if (huff)
                    ptr[planeY->some_half_array[i] * 2 + 1] = r30 | huff;
                else
                {
                    ptr[planeY->some_half_array[i] * 2 + 1] = r30;
                    mcplaneY->unk0 = decodeHuff(&state->bufTree2[0]);
                }
            }
        }
        HVQPlaneDesc *planeU = &state->planes[1];
        MCPlane *mcplaneU = &mcplanes[1];
        MCPlane *mcplaneV = &mcplanes[2];
        for (int i = 0; i < planeU->block_size_in_samples; ++i)
        {
            uint8_t *ptrU = mcplaneU->payload_ptr8;
            uint8_t *ptrV = mcplaneV->payload_ptr8;
            if (mcplaneU->unk0)
            {
                ptrU[planeU->some_half_array[i] * 2 + 1] = r30;
                ptrV[planeU->some_half_array[i] * 2 + 1] = r30;
                --mcplaneU->unk0;
            }
            else
            {
                int16_t huff = decodeHuff(&state->bufTree1[1]);
                if (huff)
                {
                    ptrU[planeU->some_half_array[i] * 2 + 1] = r30 | ((huff >> 0) & 0xF);
                    ptrV[planeU->some_half_array[i] * 2 + 1] = r30 | ((huff >> 4) & 0xF);
                }
                else
                {
                    ptrU[planeU->some_half_array[i] * 2 + 1] = r30;
                    ptrV[planeU->some_half_array[i] * 2 + 1] = r30;
                    mcplaneU->unk0 = decodeHuff(&state->bufTree2[1]);
                }
            }
        }
    }
}

static void spread_PB_descMap(SeqObj *seqobj, MCPlane mcplanes[3])
{
    uint32_t proc[2];
    uint32_t type[2];
    VideoState *state = seqobj->state;
    initMCBproc(&state->bufTree4[0], proc);
    initMCBtype(&state->bufTree4[1], type);
    for (int i = 0; i < seqobj->height; i += 8)
    {
        setMCTop(mcplanes);
        for (int j = 0; j < seqobj->width; j += 8)
        {
            getMCBtype(&state->bufTree4[1], type);
            if (type[0] == 0)
            {
                decode_PB_dc(state, mcplanes);
                decode_PB_cc(state, mcplanes, 0, type[0]);
            }
            else
            {
                reset_PB_dc(mcplanes);
                decode_PB_cc(state, mcplanes, getMCBproc(&state->bufTree4[0], proc), type[0]);
            }
            setMCNextBlk(mcplanes);
        }
        setMCDownBlk(mcplanes);
    }
}

static void resetMCHandler(VideoState *state, MCPlane mcplanes[3], void *present)
{
    for (int i = 0; i < 3; ++i)
    {
        mcplanes[i].present = present;
        mcplanes[i].payload_ptr8 = state->planes[i].payload;
        mcplanes[i].payload_ptrC = state->planes[i].payload;
        present += state->planes[i].size_in_samples;
    }
}

static void MCBlockDecDCNest(VideoState *state, MCPlane mcplanes[3])
{
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        uint8_t const *ptr = mcplanes[plane_idx].payload_ptr8;
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        uint32_t stride = plane->width_in_samples;
        int32_t line = plane->h_blocks_safe;
        for (int j = 0; j < plane->block_size_in_samples; ++j)
        {
            void *dst = mcplanes[plane_idx].top + plane->some_word_array[j];
            int32_t r30 = plane->some_half_array[j];
            uint32_t value = ptr[r30 * 2];
            uint32_t type = ptr[r30 * 2 + 1] & 0xF;
            // see also IpicBlockDec
            if (type == 0)
            {
                uint8_t top    = ptr[(r30 - line) * 2 + 1] & 0x77 ? value : ptr[(r30 - line) * 2];
                uint8_t left   = ptr[(r30 -   1) * 2 + 1] & 0x77 ? value : ptr[(r30 -   1) * 2];
                uint8_t right  = ptr[(r30 +   1) * 2 + 1] & 0x77 ? value : ptr[(r30 +   1) * 2];
                uint8_t bottom = ptr[(r30 + line) * 2 + 1] & 0x77 ? value : ptr[(r30 + line) * 2];
                WeightImBlock(dst, stride, value, top, bottom, left, right);
            }
            else if (type == 8)
            {
                dcBlock(dst, stride, value);
            }
            else
            {
                IntraAotBlock(state, dst, stride, value, type, plane_idx);
            }
        }
    }
}

static void setMCTarget(MCPlane mcplanes[3], uint32_t reference_frame)
{
    if (reference_frame == 0)
    {
        mcplanes[0].target = mcplanes[0].past;
        mcplanes[1].target = mcplanes[1].past;
        mcplanes[2].target = mcplanes[2].past;
    }
    else
    {
        mcplanes[0].target = mcplanes[0].future;
        mcplanes[1].target = mcplanes[1].future;
        mcplanes[2].target = mcplanes[2].future;
    }
}

static void getMVector(int32_t *result, BitBufferWithTree *buf, int32_t unk5)
{
    int32_t r31 = 1 << (unk5 + 5);
    int32_t r28 = decodeHuff(buf) << unk5;
    for (int i = unk5 - 1; i >= 0; --i)
        r28 += getBit(&buf->buf) << i;
    *result += r28;
    if (*result < r31)
    {
        if (*result < -r31)
            *result += r31 << 1;
    }
    else
    {
        *result -= r31 << 1;
    }
}

static void MCBlockDecMCNest(VideoState *state, MCPlane mcplanes[3], int32_t x, int32_t y)
{
    void *target = mcplanes[0].target + x/2 + (y/2 - 16)*state->planes[0].width_in_samples - 32;
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        MCPlane *mcplane = &mcplanes[plane_idx];
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        for (int i = 0; i < plane->block_size_in_samples; ++i)
        {
            uint8_t const *ptr = mcplane->payload_ptr8;
            uint8_t block_type = ptr[plane->some_half_array[i] * 2 + 1] & 0xF;
            void *dst = mcplane->top + plane->some_word_array[i];
            uint32_t stride = plane->width_in_samples;
            if (block_type == 6)
            {
                OrgBlock(state, dst, stride, plane_idx);
            }
            else
            {
                int32_t foo = x >> plane->width_shift;
                int32_t bar = y >> plane->height_shift;
                void const *src = mcplane->target + (bar >> 1) * plane->width_in_samples + (foo >> 1) + plane->some_word_array[i];
                if (block_type == 0)
                {
                    _MotionComp(dst, stride, src, stride, foo & 1, bar & 1);
                }
                else
                {
                    uint32_t strideY = state->planes[0].width_in_samples;
                    PrediAotBlock(state, dst, src, stride, block_type, target, strideY, plane_idx, foo & 1, bar & 1);
                }
            }
        }
    }
}

static void BpicPlaneDec(SeqObj *seqobj, void *present, void *past, void *future)
{
    MCPlane mcplanes[3];
    VideoState *state = seqobj->state;
    initMCHandler(state, mcplanes, present, past, future);
    spread_PB_descMap(seqobj, mcplanes);
    resetMCHandler(state, mcplanes, present);
    int32_t vec0, vec1;
    int32_t r30 = -1;
    for (int i = 0; i < seqobj->height; i += 8)
    {
        setMCTop(mcplanes);
        for (int j = 0; j < seqobj->width; j += 8)
        {
            uint8_t bits = *((uint8_t*)mcplanes[0].payload_ptr8 + 1);
            int8_t bla = (bits >> 5) & 3;
            if (bla == 0)
            {
                MCBlockDecDCNest(state, mcplanes);
            }
            else
            {
                --bla;
                if (bla != r30)
                {
                    r30 = bla;
                    setMCTarget(mcplanes, r30);
                    vec0 = 0;
                    vec1 = 0;
                }
                getMVector(&vec0, &state->bufTree3[0], state->unk6CD2[r30]);
                getMVector(&vec1, &state->bufTree3[1], state->unk6CD4[r30]);
                if (((bits >> 4) & 1) == 0)
                    MCBlockDecMCNest(state, mcplanes, j*2 + vec0, i*2 + vec1);
                else
                    MotionComp(state, mcplanes, j*2 + vec0, i*2 + vec1);
            }
            setMCNextBlk(mcplanes);
        }
        setMCDownBlk(mcplanes);
    }
}

static void HVQM4DecodeIpic(SeqObj *seqobj, uint8_t const *frame, void *present)
{
    VideoState *state = seqobj->state;
    uint8_t scale = *frame++;
    state->unk_shift = *frame++;
    frame += 2; // unused, seems to be always zero
    uint16_t nest_x = read16(frame); frame += 2;
    uint16_t nest_y = read16(frame); frame += 2;
    uint8_t const *data = frame + 0x40;
    for (int i = 0; i < 2; ++i)
    {
        setCode(&state->bufTree1[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->bufTree2[i].buf, data + read32(frame)); frame += 4;
    }
    for (int i = 0; i < 3; ++i)
    {
        setCode(&state->symbBuff[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->bufTree0[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->buf0[i],         data + read32(frame)); frame += 4;
    }
    for (int i = 0; i < 3; ++i)
    {
        setCode(&state->dc_rle[i].buf, data + read32(frame)); frame += 4;
    }
    readTree(&state->bufTree1[0], 0, 0);
    readTree(&state->bufTree2[0], 0, 0);
    readTree(&state->symbBuff[0], 1, scale);
    readTree(&state->bufTree0[0], 0, 2);

    state->boundB = +0x7F << scale;
    state->boundA = -0x80 << scale;

    Ipic_BasisNumDec(state);
    IpicDcvDec(state);
    MakeNest(state, nest_x, nest_y);

    for (int i = 0; i < 3; ++i)
    {
        IpicPlaneDec(state, i, present);
        present += state->planes[i].size_in_samples;
    }
}

static void HVQM4DecodeBpic(SeqObj *seqobj, uint8_t const *frame, void *present, void *past, void *future)
{
    VideoState *state = seqobj->state;
    state->unk_shift = frame[1];
    state->unk6CD1 = frame[0];
    state->unk6CD2[0] = frame[2];
    state->unk6CD4[0] = frame[3];
    state->unk6CD2[1] = frame[4];
    state->unk6CD4[1] = frame[5];
    // frame[6] and frame[7] are unused
    frame += 8;
    uint8_t const *data = frame + 0x44;
    for (int i = 0; i < 2; ++i)
    {
        setCode(&state->bufTree1[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->bufTree2[i].buf, data + read32(frame)); frame += 4;
    }
    for (int i = 0; i < 3; ++i)
    {
        setCode(&state->symbBuff[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->bufTree0[i].buf, data + read32(frame)); frame += 4;
        setCode(&state->buf0[i],         data + read32(frame)); frame += 4;
    }
    for (int i = 0; i < 2; ++i)
    {
        setCode(&state->bufTree3[i].buf, data + read32(frame)); frame += 4;
    }
    for (int i = 0; i < 2; ++i)
    {
        setCode(&state->bufTree4[1 - i].buf, data + read32(frame)); frame += 4;
    }
    readTree(&state->bufTree1[0], 0, 0);
    readTree(&state->bufTree2[0], 0, 0);
    readTree(&state->symbBuff[0], 1, state->unk6CD1);
    readTree(&state->bufTree0[0], 0, 2);
    readTree(&state->bufTree3[0], 1, 0);
    readTree(&state->bufTree4[1], 0, 0);

    state->boundB = +0x7F << state->unk6CD1;
    state->boundA = -0x80 << state->unk6CD1;

    BpicPlaneDec(seqobj, present, past, future);
}

static void HVQM4DecodePpic(SeqObj *seqobj, uint8_t const *frame, void *present, void *past)
{
    HVQM4DecodeBpic(seqobj, frame, present, past, present);
}

// I/P frames shift the past/future reference window
// B frames are never used as references
enum FrameType
{
    I_FRAME = 0x10,
    P_FRAME = 0x20,
    B_FRAME = 0x30,
};

#ifndef NATIVE
#define YOLO_INCLUDE
#include "yoloader.c"
#endif


static void decode_video(Player *player, FILE *infile, uint32_t gop_start, uint16_t frame_type, uint32_t frame_size)
{
    // getBit() and getByte() overread by up to 3 bytes
    uint32_t overread = 3;
    uint8_t *frame = malloc(frame_size + overread);
    fread(frame, frame_size, 1, infile);

    uint32_t disp_id = read32(frame);

    // swap past and future
    if (frame_type != B_FRAME)
    {
        void *tmp = player->past;
        player->past = player->future;
        player->future = tmp;
    }

    SeqObj *seqobj = &player->seqobj;
    switch (frame_type)
    {
    case I_FRAME: HVQM4DecodeIpic(seqobj, frame + 4, player->present);                               break;
    case P_FRAME: HVQM4DecodePpic(seqobj, frame + 4, player->present, player->past);                 break;
    case B_FRAME: HVQM4DecodeBpic(seqobj, frame + 4, player->present, player->past, player->future); break;
    default:
        fprintf(stderr, "unknown video frame type 0x%x\n", frame_type);
        exit(EXIT_FAILURE);
    }
    free(frame);

    //if (frame_type == I_FRAME)
    {
        char name[50];
        sprintf(name, "output/video_rgb_%u.ppm", gop_start + disp_id);
        //char type = "ipb"[(frame_type >> 4) - 1];
        //sprintf(name, "output/video_rgb_%u_%u_%c.ppm", gop, disp_id, type);
        //printf("writing frame to %s...\n", name);
        dumpRGB(player, name);
        //sprintf(name, "output/video_yuv_%u_%u_%c.ppm", gop, disp_id, type);
        //dumpYUV(player, name);
    }

    // swap present and future
    if (frame_type != B_FRAME)
    {
        void *tmp = player->present;
        player->present = player->future;
        player->future = tmp;
    }
}

/* stream structure */

const char HVQM4_13_magic[16] = "HVQM4 1.3";
const char HVQM4_15_magic[16] = "HVQM4 1.5";

typedef struct
{
    enum
    {
        HVQM4_13,
        HVQM4_15,
    } version;

    uint32_t header_size;   /* 0x10-0x13 */
    uint32_t body_size;     /* 0x14-0x17 */
    uint32_t blocks;        /* 0x18-0x1B */
    uint32_t video_frames;  /* 0x1C-0x1F */
    uint32_t audio_frames;  /* 0x20-0x23 */
    uint32_t usec_per_frame;/* 0x24-0x27 (33366, 33367, 40000) */
    uint32_t max_frame_sz;  /* 0x28-0x2B */
    uint32_t unk2C;         /* 0x2C-0x2F (0) */
    uint32_t audio_frame_sz;/* 0x30-0x33 */
    uint16_t hres;          /* 0x34-0x35 */
    uint16_t vres;          /* 0x36-0x37 */
    uint8_t  hsamp;         /* 0x38 */
    uint8_t  vsamp;         /* 0x39 */
    uint8_t  video_mode;    /* 0x3A (0 or 0x12) */
    uint8_t  unk3B;         /* 0x3B (0) */
    uint8_t  audio_channels;/* 0x3C */
    uint8_t  audio_bitdepth;/* 0x3D */
    uint16_t pad;           /* 0x3E-0x3F (0) */
    uint32_t audio_srate;   /* 0x40-0x43 */
} HVQM4_header;

static void load_header(HVQM4_header *header, uint8_t *raw_header)
{
    /* check MAGIC */
    if (!memcmp(HVQM4_13_magic, &raw_header[0], 16))
    {
        header->version = HVQM4_13;
    }
    else if (!memcmp(HVQM4_15_magic, &raw_header[0], 16))
    {
        header->version = HVQM4_15;
    }
    else
    {
        fprintf(stderr, "does not appear to be a HVQM4 file\n");
        exit(EXIT_FAILURE);
    }

    header->header_size = read32(&raw_header[0x10]);
    header->body_size = read32(&raw_header[0x14]);
    header->blocks = read32(&raw_header[0x18]);
    header->video_frames = read32(&raw_header[0x1C]);
    header->audio_frames = read32(&raw_header[0x20]);
    header->usec_per_frame = read32(&raw_header[0x24]);
    header->max_frame_sz = read32(&raw_header[0x28]);
    header->unk2C = read32(&raw_header[0x2C]);
    header->audio_frame_sz = read32(&raw_header[0x30]);
    header->hres = read16(&raw_header[0x34]);
    header->vres = read16(&raw_header[0x36]);
    header->hsamp = raw_header[0x38];
    header->vsamp = raw_header[0x39];
    header->video_mode = raw_header[0x3A];
    header->unk3B = raw_header[0x3B];
    header->audio_channels = raw_header[0x3C];
    header->audio_bitdepth = raw_header[0x3D];
    header->pad = read16(&raw_header[0x3E]);
    header->audio_srate = read32(&raw_header[0x40]);

    expect32_text(0x44, header->header_size, "header size");
    /* no check for body size yet */
    if (header->blocks == 0)
    {
        fprintf(stderr, "zero blocks\n");
        exit(EXIT_FAILURE);
    }
#if 0
    if (header->audio_frames != 0)
    {
        if (header->audio_srate == 0 || header->audio_frame_sz == 0 || header->audio_channels == 0)
        {
            fprintf(stderr, "expected nonzero audio srate and frame size\n");
            exit(EXIT_FAILURE);
        }
    }
#endif
    /* no check for video frame count */
    /*
    if (header->unk24 != 0x8257 && header->unk24 != 0x8256)
    {
        expect32_imm(0x8257, header->unk24, 0x24);
    }
    */
    expect32_imm(0, header->unk2C, 0x2C);
    //expect32_text(0x650, header->audio_frame_sz, "audio frame size");
    /* no check for hres, vres */
    if (header->video_mode != 0 && header->video_mode != 0x12)
    {
        expect32_imm(0, header->video_mode, 0x3A);
    }
    expect32_imm(0, header->unk3B, 0x3B);
    //expect32_imm(0x02100000, header->unk3C, 0x3C);
    /* no check for srate, can be 0 */
}

static void display_header(HVQM4_header *header)
{
    switch (header->version)
    {
        case HVQM4_13:
            printf("HVQM4 1.3\n");
            break;
        case HVQM4_15:
            printf("HVQM4 1.5\n");
            break;
    }
    printf("Header size: 0x%"PRIx32"\n", header->header_size);
    printf("Body size:   0x%"PRIx32"\n", header->body_size);
    printf("Max frame:   0x%"PRIx32"\n", header->max_frame_sz);
    printf("Resolution:  %"PRIu32" x %"PRIu32"\n", header->hres, header->vres);
    printf("Chroma subsampling: %u x %u\n", header->hsamp, header->vsamp);
    printf("Video mode: %u\n", header->video_mode);
    printf("µs/frame:    %"PRIu32"\n", header->usec_per_frame);
    printf("%d GOPs\n", header->blocks);
    printf("%d Video frames\n", header->video_frames);
    printf("%d Audio frames\n", header->audio_frames);
    printf("Sample rate: %"PRIu32" Hz\n", header->audio_srate);
    printf("Audio frame size: 0x%"PRIx32"\n", header->audio_frame_sz);
    printf("Audio channels: %u\n", header->audio_channels);
    printf("\n");
}

static void put_32bitLE(uint8_t * buf, uint32_t v)
{
    for (unsigned int i = 0; i < 4; i++)
    {
        buf[i] = v & 0xFF;
        v >>= 8;
    }
}

static void put_16bitLE(uint8_t * buf, uint16_t v)
{
    for (unsigned int i = 0; i < 2; i++)
    {
        buf[i] = v & 0xFF;
        v >>= 8;
    }
}

/* make a header for 16-bit PCM .wav */
/* buffer must be 0x2c bytes */
static void make_wav_header(uint8_t * buf, int32_t sample_count, int32_t sample_rate, int channels) {
    size_t bytecount;

    bytecount = sample_count*channels*2;

    /* RIFF header */
    memcpy(buf+0, "RIFF", 4);
    /* size of RIFF */
    put_32bitLE(buf+4, (int32_t)(bytecount+0x2c-8));

    /* WAVE header */
    memcpy(buf+8, "WAVE", 4);

    /* WAVE fmt chunk */
    memcpy(buf+0xc, "fmt ", 4);
    /* size of WAVE fmt chunk */
    put_32bitLE(buf+0x10, 0x10);

    /* compression code 1=PCM */
    put_16bitLE(buf+0x14, 1);

    /* channel count */
    put_16bitLE(buf+0x16, channels);

    /* sample rate */
    put_32bitLE(buf+0x18, sample_rate);

    /* bytes per second */
    put_32bitLE(buf+0x1c, sample_rate*channels*2);

    /* block align */
    put_16bitLE(buf+0x20, (int16_t)(channels*2));

    /* significant bits per sample */
    put_16bitLE(buf+0x22, 2*8);

    /* PCM has no extra format bytes, so we don't even need to specify a count */

    /* WAVE data chunk */
    memcpy(buf+0x24, "data", 4);
    /* size of WAVE data chunk */
    put_32bitLE(buf+0x28, (int32_t)bytecount);
}

static void decv_init(Player *player)
{
    // HACK
    uint32_t res = player->seqobj.width * player->seqobj.height;
    uint32_t sample_size = player->seqobj.h_samp * player->seqobj.v_samp;
    uint32_t picsize = (res * (sample_size + 2)) / sample_size;

    player->past    = malloc(picsize);
    player->present = malloc(picsize);
    player->future  = malloc(picsize);
}

#ifndef HVQM4_NOMAIN
int main(int argc, char **argv)
{
    printf("h4m 'HVQM4 1.3/1.5' decoder 0.4 by flacs/hcs\n\n");
    if (argc != 3)
    {
        fprintf(stderr, "usage: %s file.h4m output.wav\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* open input */
    FILE *infile = fopen(argv[1], "rb");
    if (!infile)
    {
        fprintf(stderr, "failed opening %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }

    /* read in header */
    uint8_t raw_header[0x44];
    if (0x44 != fread(&raw_header, 1, 0x44, infile))
    {
        fprintf(stderr, "failed reading header");
        exit(EXIT_FAILURE);
    }

    /* load up and check header */
    HVQM4_header header;
    load_header(&header, raw_header);
    display_header(&header);

    /* open output */
    FILE *outfile = fopen(argv[2], "wb");
    if (!outfile)
    {
        fprintf(stderr, "error opening %s\n", argv[2]);
        exit(EXIT_FAILURE);
    }

    /* fill in the space that we'll put the header in later */
    uint8_t riff_header[0x2c];
    memset(riff_header, 0, sizeof(riff_header));
    if (0x2c != fwrite(riff_header, 1, 0x2c, outfile))
    {
        fprintf(stderr, "error writing riff header\n");
        exit(EXIT_FAILURE);
    }

#ifndef NATIVE
    load_library("RM_DLL.elf");
    pHVQM4InitDecoder();
#endif

    HVQM4InitDecoder();

    Player player;
#ifdef NATIVE
    HVQM4InitSeqObj(&player.seqobj, (VideoInfo*)&header.hres);
    VideoState *state = malloc(HVQM4BuffSize(&player.seqobj));
    HVQM4SetBuffer(&player.seqobj, state);
#else
    pHVQM4InitSeqObj(&player.seqobj, (VideoInfo*)&header.hres);
    VideoState *state = malloc(pHVQM4BuffSize(&player.seqobj));
    pHVQM4SetBuffer(&player.seqobj, state);
#endif
    decv_init(&player);

    /* parse blocks */
    uint32_t block_count = 0;
    uint32_t total_aud_frames = 0;
    uint32_t total_vid_frames = 0;
    uint32_t last_block_start = ftell(infile);
    uint32_t total_sample_count = 0;
    while (block_count < header.blocks)
    {
        const long block_start = ftell(infile);
        //expect32(ftell(infile) - last_block_start, infile);
        get32(infile);
        last_block_start = block_start;
        const uint32_t expected_block_size = get32(infile);
        const uint32_t expected_vid_frame_count = get32(infile);
        const uint32_t expected_aud_frame_count = get32(infile);
        expect32(0x01000000, infile);   /* EOS marker? */
        //get32(infile);   /* EOS marker? */
        const long data_start = ftell(infile);

        block_count ++;
#ifdef VERBOSE_PRINT
        printf("\n\nblock %u/%u starts at 0x%lx, length 0x%"PRIx32"\n", (int)block_count, header.blocks, block_start, expected_block_size);
#endif

        /* parse frames */
        struct audio_state audio_state;
        int first_aud=1;
        uint32_t vid_frame_count = 0, aud_frame_count = 0;
        int block_sample_count =0;

        audio_state.ch = calloc(header.audio_channels, sizeof(*audio_state.ch));
        while (aud_frame_count < expected_aud_frame_count ||
               vid_frame_count < expected_vid_frame_count)
        {
            const uint16_t frame_id1 = get16(infile);
            const uint16_t frame_id2 = get16(infile);
            const uint32_t frame_size = get32(infile);

#ifdef VERBOSE_PRINT
            printf("\nframe id 0x%04"PRIx16",0x%04"PRIx16" ",frame_id1,frame_id2);
            printf("size 0x%08"PRIx32"\n", frame_size);
#endif

            if (frame_id1 == 1)
            {
                /* video */
                vid_frame_count ++;
#ifdef VERBOSE_PRINT
                printf("video frame %d/%d\n", (int)vid_frame_count, (int)expected_vid_frame_count);
#endif
                decode_video(&player, infile, total_vid_frames, frame_id2, frame_size);
            }
            else if (frame_id1 == 0)
            {
                /* audio */
                aud_frame_count ++;
                total_aud_frames ++;
                const long audio_started = ftell(infile);
                const uint32_t samples = get32(infile);
                block_sample_count += samples;
#ifdef VERBOSE_PRINT
                printf("0x%lx: audio frame %d/%d (%d) (%d samples)\n", (unsigned long)audio_started, (int)aud_frame_count, (int)expected_aud_frame_count, (int)total_aud_frames, samples);
#endif
#if 1
                decode_audio(&audio_state, first_aud, samples, infile, outfile, header.audio_channels);
                first_aud = 0;
                long bytes_done_unto = ftell(infile) - audio_started;
                if (bytes_done_unto > frame_size)
                {
                    fprintf(stderr, "processed 0x%lx bytes, should have done 0x%"PRIx32"\n",
                        bytes_done_unto, frame_size);
                    exit(EXIT_FAILURE);
                }
                else if (bytes_done_unto < frame_size)
                {
                    while (bytes_done_unto < frame_size)
                    {
                        get8(infile);
                        bytes_done_unto ++;
                    }
                }
#else
                printf("skipping\n");
                seek_past(frame_size - 4, infile);
#endif
            }
            else
            {
                fprintf(stderr, "unexpected frame id 0x%04X 0x%04X at %08lx\n", frame_id1, frame_id2, (unsigned long)(ftell(infile)-8));
                exit(EXIT_FAILURE);
                seek_past(frame_size, infile);
            }
        }
        free(audio_state.ch);

        if (expected_aud_frame_count != aud_frame_count ||
            expected_vid_frame_count != vid_frame_count)
        {
            fprintf(stderr, "frame count mismatch\n");
            exit(EXIT_FAILURE);
        }

        total_vid_frames += vid_frame_count;

#ifdef VERBOSE_PRINT
        printf("block %d ended at 0x%lx (%d samples)\n", (int)block_count, ftell(infile), block_sample_count);
#endif
        if (ftell(infile) != (data_start+expected_block_size))
        {
            fprintf(stderr, "block size mismatch\n");
            exit(EXIT_FAILURE);
        }
        total_sample_count += block_sample_count;
    }

    free(player.seqobj.state);
    free(player.past);
    free(player.present);
    free(player.future);

    if (total_aud_frames != header.audio_frames ||
        total_vid_frames != header.video_frames)
    {
        fprintf(stderr, "total frame count mismatch\n");
        exit(EXIT_FAILURE);
    }

    printf("%"PRIu32" samples\n", total_sample_count);

    // generate header
    make_wav_header(riff_header, total_sample_count, header.audio_srate, header.audio_channels);
    fseek(outfile, 0, SEEK_SET);
    if (0x2c != fwrite(riff_header, 1, 0x2c, outfile))
    {
        fprintf(stderr, "error rewriting riff header\n");
        exit(EXIT_FAILURE);
    }
    if (EOF == fclose(outfile))
    {
        fprintf(stderr, "error finishing output\n");
        exit(EXIT_FAILURE);
    }

    printf("Done!\n");
}
#endif
