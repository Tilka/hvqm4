#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

/* .h4m (HVQM4 1.3/1.5) audio decoder 0.3 by hcs */

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

static void expect32_text(uint32_t expected, uint32_t actual, char *name)
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

    if (sample_count != fwrite(samples, sizeof(int16_t)*channels, sample_count, outfile))
    {
        fprintf(stderr, "error writing output\n");
        exit(EXIT_FAILURE);
    }

    free(samples);
}

/* video stuff */

static uint8_t clipTable[0x200];
static int32_t divTable[0x10];
static int32_t mcdivTable[0x200];

static void init_global_constants()
{
    for (int i = 0, j = -0x80; i < 0x200; ++i, ++j)
        clipTable[i] = j < 0 ? 0 : (j > 0xff ? 0xff : (j & 0xff));
    divTable[0] = 0;
    mcdivTable[0] = 0;
    for (int i = 1; i < 0x10; ++i)
        divTable[i] = 0x1000 / (i * 16) * 16;
    for (int i = 1; i < 0x200; ++i)
        mcdivTable[i] = 0x1000 / i;
}

void HVQM4InitDecoder()
{
    init_global_constants();
}

// done
static void dcBlock(uint8_t *dst, uint32_t stride, uint8_t value)
{
    printf("  dcBlock()\n");
    for (int i = 0; i < 4; ++i)
    {
        dst[0] = value;
        dst[1] = value;
        dst[2] = value;
        dst[3] = value;
        dst += stride;
    }
}

static void WeightImBlock(uint8_t *dst, uint32_t stride, uint8_t unk5, uint8_t unk6, uint8_t unk7, uint8_t unk8, uint8_t unk9)
{
    //printf("  WeightImBlock(dst=%p, stride=%u, %u, %u,%u,%u,%u): ", dst, stride, unk5, unk6, unk7, unk8, unk9);
    int32_t diff0 = unk6 - unk7;
    int32_t diff1 = unk8 - unk9;
    int32_t r29 = diff0 + diff1;
    int32_t r28 = diff0 - diff1;

    int32_t r27 = unk5 << 1;
    int32_t r31 = (unk5 << 3) + 4;

    int32_t r24 = (unk6 + unk8) - r27;
    int32_t r23 = (unk6 + unk9) - r27;
    int32_t r22 = (unk7 + unk9) - r27;
    int32_t r21 = (unk8 + unk7) - r27;

    int32_t r20 = unk6 - unk8;
    int32_t r19 = unk6 - unk9;
    int32_t r18 = unk7 - unk9;
    int32_t r17 = unk7 - unk8;

    dst[0] = clipTable[((r29 + r24 + r31) >> 3) + 0x80];
    dst[1] = clipTable[((r29 + r20 + r31) >> 3) + 0x80];
    dst[2] = clipTable[((r28 + r19 + r31) >> 3) + 0x80];
    dst[3] = clipTable[((r28 + r23 + r31) >> 3) + 0x80];
    //printf("%u %u %u %u  ", dst[0], dst[1], dst[2], dst[3]);

    dst += stride;

    dst[0] = clipTable[((r31 + r29 - r20) >> 3) + 0x80];
    dst[1] = clipTable[((r31 - r22      ) >> 3) + 0x80];
    dst[2] = clipTable[((r31 - r21      ) >> 3) + 0x80];
    dst[3] = clipTable[((r31 + r28 - r19) >> 3) + 0x80];
    //printf("%u %u %u %u  ", dst[0], dst[1], dst[2], dst[3]);

    dst += stride;

    dst[0] = clipTable[((r31 - r28 - r17) >> 3) + 0x80];
    dst[1] = clipTable[((r31 - r23      ) >> 3) + 0x80];
    dst[2] = clipTable[((r31 - r24      ) >> 3) + 0x80];
    dst[3] = clipTable[((r31 - r29 - r18) >> 3) + 0x80];
    //printf("%u %u %u %u  ", dst[0], dst[1], dst[2], dst[3]);

    dst += stride;

    dst[0] = clipTable[((r31 - r28 + r21) >> 3) + 0x80];
    dst[1] = clipTable[((r31 - r28 + r17) >> 3) + 0x80];
    dst[2] = clipTable[((r31 - r29 + r18) >> 3) + 0x80];
    dst[3] = clipTable[((r31 - r29 + r22) >> 3) + 0x80];
    //printf("%u %u %u %u\n", dst[0], dst[1], dst[2], dst[3]);
}

typedef struct
{
    uint32_t pos;
    int32_t tree_unk4;
    uint32_t array0[0x200]; // 8+
    uint32_t array1[0x200]; // 0x808+
} Tree;

typedef struct
{
    void const *ptr;   // 0-3
    uint32_t buf_unk4; // 4-7
    uint32_t value;    // 8-B
    int32_t bit;       // C-F
} BitBuffer;

typedef struct
{
    BitBuffer buf; // 0-F
    Tree *tree;    // 0x10-0x13
} BitBufferWithTree;

typedef struct
{
    // size: 0x38 (56)
    // lots of little things in here
    void *ptr0; // 0-3 beginning of the plane including the border
    void *ptr4; // 4-7 beginning of the non-border plane
    uint16_t h_blocks; // 8-9
    uint16_t v_blocks; // A-B
    uint16_t h_blocks_safe; // C-D
    uint16_t v_blocks_safe; // E-F
    // [...]
    uint16_t width_in_samples; // 28-29
    uint16_t height_in_samples; // 2A-2B
    uint32_t size_in_samples; // 2C-2F
    uint8_t width_shift; // 30
    uint8_t height_shift; // 31
    uint8_t h_samp_per_block; // 32
    uint8_t v_samp_per_block; // 33
    uint8_t block_size_in_samples; // 34
} HVQPlaneDesc;

typedef struct
{
    uint8_t yuvbuf[640*480*3]; // HACK
    // size: 0x6CD8 (plane data comes right after)
    HVQPlaneDesc planes[3]; // 0x00 - 0xA8
    Tree trees[6];
    BitBufferWithTree symbBuff[3];
    BitBufferWithTree huffBuff[3];
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
    uint32_t boundA; // 0x6CCC
    uint32_t boundB; // 0x6CC8
    uint8_t unk_shift; // 0x6CD0
    uint8_t maybe_padding[7]; // 0x6CD1-0x6CD7
} VideoState;

typedef struct
{
    VideoState *state;
    uint16_t width;
    uint16_t height;
    uint8_t h_samp;
    uint8_t v_samp;
} SeqObj;

typedef struct
{
    uint16_t hres;
    uint16_t vres;
    uint8_t h_samp;
    uint8_t v_samp;
    uint8_t video_mode;
} VideoInfo;

static void OrgBlock(VideoState *state, uint8_t *dst, uint32_t stride, uint32_t bufIndex)
{
    printf("OrgBlock\n");
    void const **p = &state->buf0[bufIndex].ptr;
    for (int i = 0; i < 4; ++i)
    {
        *(uint32_t*)(dst + i * stride) = read32(*p + i);
        (*p) += 4;
    }
}

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

static int16_t getByte(BitBuffer *buf)
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
        int16_t byte = getByte(src);
        int16_t value = byte;
        if (readTree_signed && value > 0x7F)
            value -= 0x100;
        value <<= readTree_scale;
        dst->array0[byte] = (uint32_t)value;
        return byte;
    }
    else
    {
        uint32_t pos = dst->pos++;
        dst->array0[pos] = (uint32_t)_readTree(dst, src);
        dst->array1[pos] = (uint32_t)_readTree(dst, src);
        return (int16_t)pos;
    }
}

static void readTree(BitBufferWithTree *buf, uint32_t isSigned, uint32_t scale)
{
    readTree_signed = isSigned;
    readTree_scale = scale;
    Tree *tree = buf->tree;
    tree->pos = 0x100;
    if (buf->buf.buf_unk4 == 0)
        tree->tree_unk4 = 0;
    else
        tree->tree_unk4 = _readTree(tree, &buf->buf);
}

static uint32_t decodeHuff(BitBufferWithTree *buf)
{
    Tree *tree = buf->tree;
    int32_t tree_unk4 = tree->tree_unk4;
    while (tree_unk4 >= 0x100)
        tree_unk4 = tree->array0[(getBit(&buf->buf) << 9) + tree_unk4];
    return tree->array0[tree_unk4];
}

static int32_t decodeSOvfSym(BitBufferWithTree *buf, int32_t a, int32_t b)
{
    int32_t sum = 0;
    int32_t value;
    do
    {
        value = decodeHuff(buf);
        sum += value;
    } while (a >= value && value >= b);
    return sum;
}

static uint32_t GetAotBasis(VideoState *state, uint8_t dst[4][4], int32_t *sum, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
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
    //printf("bits=0x%04X, step=0x%02X, stride=0x%02X, nest_stride=0x%02X\n", bits, step, stride, nest_stride);
    uint8_t min = *nest_data;
    uint8_t max = *nest_data;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            uint8_t nest_value = nest_data[i * stride + j * step];
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

static int32_t GetAot1(VideoState *state, int32_t result[4][4], uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    uint8_t byte_result[4][4];
    int32_t dummy = 0;
    uint32_t factor = GetAotBasis(state, byte_result, &dummy, nest_data, nest_stride, plane_idx);
    int32_t sum = 0;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            sum += result[i][j] = factor * byte_result[i][j];
    return sum >> 4;
}

static int32_t GetAotSum(VideoState *state, int32_t result[4][4], uint8_t some_counter, uint8_t const *nest_data, uint32_t nest_stride, uint32_t plane_idx)
{
    // TODO
    while (some_counter--)
    {
        // GetAotBasis(state, );
    }
    return 0;
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

static void setHVQPlaneDesc(SeqObj *seqobj, uint8_t plane_idx, uint8_t h_samp, uint8_t v_samp)
{
    HVQPlaneDesc *plane = &seqobj->state->planes[plane_idx];
    plane->width_shift = h_samp == 2 ? 1 : 0;
    plane->width_in_samples = seqobj->width >> plane->width_shift;
    plane->height_shift = v_samp == 2 ? 1 : 0;
    plane->height_in_samples = seqobj->height >> plane->height_shift;
    plane->size_in_samples = plane->width_in_samples * plane->height_in_samples;
    plane->h_samp_per_block = 2 >> plane->width_shift;
    plane->v_samp_per_block = 2 >> plane->height_shift;
    plane->block_size_in_samples = plane->h_samp_per_block * plane->v_samp_per_block;
    plane->h_blocks = seqobj->width / (h_samp * 4);
    plane->v_blocks = seqobj->height / (v_samp * 4);
    plane->h_blocks_safe = plane->h_blocks + 2;
    plane->v_blocks_safe = plane->v_blocks + 2;
    // TODO: a bit more stuff i don't need yet
}

// HACK: assumes 4:2:0
static void dumpYUV(SeqObj *seqobj, char const *path)
{
    FILE *f = fopen(path, "wb+");
    fprintf(f, "P2\n%u %u\n255\n", seqobj->width, seqobj->height * 2);
    uint8_t const *p = seqobj->state->yuvbuf;
    printf("dumping memory from %p\n", p);
    for (int plane = 0; plane < 2; ++plane)
    {
        for (int i = 0; i < seqobj->height; ++i)
        {
            for (int j = 0; j < seqobj->width; ++j)
            {
                if (plane == 0 || j < seqobj->width / 2)
                    fprintf(f, "%u ", *p++);
                else
                    fprintf(f, "0 ");
            }
            fputs("\n", f);
        }
    }
    fclose(f);
}

static void dumpPlanes(VideoState *state, char const *prefix)
{
    printf("dumping %s\n", prefix);
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        char path[128];
        snprintf(path, 128, "%s_%c.ppm", prefix, "yuv"[plane_idx]);
        FILE *f = fopen(path, "wb+");
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        fprintf(f, "P2\n%u %u\n255\n", plane->h_blocks_safe, plane->v_blocks_safe);
        uint8_t const *p = plane->ptr0;
        for (int i = 0; i < plane->v_blocks_safe; ++i)
        {
            for (int j = 0; j < plane->h_blocks_safe; ++j)
            {
                fprintf(f, "%u ", *p++);
                //fprintf(f, "%u ", *p++);
                ++p;
            }
            fputs("\n", f);
        }
        fclose(f);
    }
}

static void set_border(uint8_t *dst)
{
    dst[0] = 0x7F;
    dst[1] = 0xFF;
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

    state->huffBuff[0].tree = &state->trees[1]; // reuse!
    state->huffBuff[1].tree = &state->trees[1]; //
    state->huffBuff[2].tree = &state->trees[1]; //

    state->bufTree0[0].tree = &state->trees[2];
    state->bufTree0[1].tree = &state->trees[2];
    state->bufTree0[2].tree = &state->trees[2];

    state->bufTree3[0].tree = &state->trees[4];
    state->bufTree3[1].tree = &state->trees[4];

    state->bufTree4[0].tree = &state->trees[5];
    state->bufTree4[1].tree = &state->trees[5];

    void *plane_data = workbuff + sizeof(VideoState);
    for (int i = 0; i < 3; ++i)
    {
        HVQPlaneDesc *plane = &state->planes[i];
        plane->ptr0 = plane_data;
        plane->ptr4 = plane_data + plane->h_blocks_safe * sizeof(uint16_t) + 2;
        plane_data += plane->h_blocks_safe * plane->v_blocks_safe * sizeof(uint16_t);

        // set horizontal borders
        void *ptr = plane->ptr0;
        for (uint32_t i = plane->h_blocks_safe; i; --i)
        {
            set_border(ptr);
            ptr += sizeof(uint16_t);
        }

        ptr = plane_data;
        for (uint32_t i = plane->h_blocks_safe; i; --i)
        {
            ptr -= sizeof(uint16_t);
            set_border(ptr);
        }

        // set vertical borders
        ptr = plane->ptr0 + plane->h_blocks_safe * sizeof(uint16_t);
        for (uint32_t i = plane->v_blocks_safe - 2; i; --i)
        {
            set_border(ptr);
            ptr += plane->h_blocks_safe * sizeof(uint16_t);
        }

        ptr = plane->ptr0 + plane->h_blocks_safe * sizeof(uint16_t) * 2 - sizeof(uint16_t);
        for (uint32_t i = plane->v_blocks_safe - 2; i; --i)
        {
            set_border(ptr);
            ptr += plane->h_blocks_safe * sizeof(uint16_t);
        }
    }
}

static uint32_t getDeltaDC(VideoState *state, uint32_t plane_idx, uint32_t *ptr)
{
    //printf("getDeltaDC(plane=%u, *ptr=%u): ", plane_idx, *ptr);
    if (*ptr == 0)
    {
        uint32_t symbol = decodeSOvfSym(&state->symbBuff[plane_idx], state->boundA, state->boundB);
        if (symbol == 0)
            *ptr = decodeHuff(&state->huffBuff[plane_idx]);
        return symbol;
    }
    else
    {
        --(*ptr);
        //puts("");
        return 0;
    }
}

static void setCode(BitBuffer *dst, void const *src)
{
    dst->buf_unk4 = read32(src);
    dst->ptr = dst->buf_unk4 ? src + 4 : NULL;
    dst->bit = -1;
}

static void Ipic_BasisNumDec(VideoState *state)
{
    uint8_t *y_dst = state->planes[0].ptr4;
    uint32_t y_v_blocks = state->planes[0].v_blocks;
    BitBufferWithTree *y_tree1 = &state->bufTree1[0];
    BitBufferWithTree *y_tree2 = &state->bufTree2[0];
    uint32_t value = 0;
    while (y_v_blocks--)
    {
        uint32_t y_h_blocks = state->planes[0].h_blocks;
        while (y_h_blocks--)
        {
            if (value)
            {
                y_dst[1] = 0;
                y_dst += 2;
                --value;
            }
            else
            {
                int16_t x = decodeHuff(y_tree1) & 0xffff;
                if (x == 0)
                    value = decodeHuff(y_tree2);
                y_dst[1] = x & 0xff;
                y_dst += 2;
            }
        }
        y_dst += 4;
    }

    uint8_t *u_dst = state->planes[1].ptr4;
    uint8_t *v_dst = state->planes[2].ptr4;
    uint32_t uv_v_blocks = state->planes[1].v_blocks;
    BitBufferWithTree *uv_tree1 = &state->bufTree1[1];
    BitBufferWithTree *uv_tree2 = &state->bufTree2[1];
    value = 0;
    while (uv_v_blocks--)
    {
        uint32_t uv_h_blocks = state->planes[1].h_blocks;
        while (uv_h_blocks--)
        {
            if (value)
            {
                u_dst[1] = 0;
                v_dst[1] = 0;
                u_dst += 2;
                v_dst += 2;
                --value;
            }
            else
            {
                int16_t x = decodeHuff(uv_tree1) & 0xffff;
                if (x == 0)
                    value = decodeHuff(uv_tree2);
                u_dst[1] = (x >> 0) & 0xf;
                v_dst[1] = (x >> 4) & 0xf;
                u_dst += 2;
                v_dst += 2;
            }
        }
        u_dst += 4;
        v_dst += 4;
    }
}

static void IpicDcvDec(VideoState *state)
{
    for (int plane_idx = 0; plane_idx < 3; ++plane_idx)
    {
        HVQPlaneDesc *plane = &state->planes[plane_idx];
        uint32_t x = 0;
        uint32_t v_blocks = plane->v_blocks;
        uint8_t *ptr = plane->ptr4;
        while (v_blocks--)
        {
            uint8_t const *ptr2 = ptr - plane->h_blocks_safe * 2;
            uint8_t y = ptr2[0];
            uint32_t h_blocks = plane->h_blocks;
            while (h_blocks--)
            {
                y += getDeltaDC(state, plane_idx, &x);
                //printf("y=%u\n", y);
                *ptr = y;
                y = (ptr2[2] + y + 1) >> 1;
                ptr += 2;
                ptr2 += 2;
            }
            // skip adjacent vertical borders
            ptr += 2 * sizeof(uint16_t);
        }
        dumpPlanes(state, "filled");
    }
}

static void MakeNest(VideoState *state, uint16_t unk4, uint16_t unk5)
{
    HVQPlaneDesc *y_plane = &state->planes[0];
    uint8_t *ptr = y_plane->ptr4 + (y_plane->h_blocks_safe * unk5 + unk4) * 2;

    uint32_t r19, r20;
    int32_t r23, r24, r25, r26;

    if (y_plane->h_blocks < state->h_nest_size)
    {
        r24 = y_plane->h_blocks;
        r26 = state->h_nest_size - y_plane->h_blocks;
        if (r26 > y_plane->h_blocks)
            r26 = y_plane->h_blocks;
        r20 = state->h_nest_size - (r24 + r26);
    }
    else
    {
        r24 = state->h_nest_size;
        r20 = 0;
        r26 = 0;
    }

    if (y_plane->v_blocks < state->v_nest_size)
    {
        r23 = y_plane->v_blocks;
        r25 = state->v_nest_size - y_plane->v_blocks;
        if (r25 > y_plane->v_blocks)
            r25 = y_plane->v_blocks;
        r19 = state->v_nest_size - (r23 + r25);
    }
    else
    {
        r23 = state->v_nest_size;
        r19 = 0;
        r25 = 0;
    }

    uint8_t *nest = state->nest_data;
    for (int i = 0; i < r23; ++i)
    {
        uint8_t const *p = ptr;
        for (int j = 0; j < r24; ++j)
        {
            *nest++ = (p[0] >> 4) & 0xF;
            p += 2;
        }
        for (int j = 0; j < r26; ++j)
        {
            p -= 2;
            *nest++ = (p[0] >> 4) & 0xF;
        }
        for (int j = 0; j < r20; ++j)
            *nest++ = 0;
        ptr += y_plane->h_blocks_safe * 2;
    }

    uint8_t const *nest2 = nest - state->h_nest_size;
    for (int i = 0; i < r25; ++i)
    {
        for (int j = 0; j < state->h_nest_size; ++j)
            *nest++ = nest2[j];
        nest2 -= state->h_nest_size;
    }

    for (int i = 0; i < r19; ++i)
    {
        for (int j = 0; j < state->h_nest_size; ++j)
            *nest++ = 0;
    }
}

typedef struct
{
    uint32_t plane_idx;
    void *ptr4;
    void *ptr8;
    void *ptrC;
    uint8_t unk10;
    uint8_t unk11;
    uint8_t unk12;
    uint8_t unk13;
    uint8_t unk14;
} StackState;

static void IntraAotBlock(VideoState *state, uint8_t *dst, uint32_t stride, uint8_t unk6, uint8_t block_type, uint32_t plane_idx)
{
    printf("  IntraAotBlock(unk6=%u, block_type=%u)\n", unk6, block_type);
    if (block_type == 6)
    {
        OrgBlock(state, dst, stride, plane_idx);
        return;
    }
#if 1
    for (int i = 0; i < 70*38; ++i)
    {
        printf("%02X ", state->nest_data[i]);
        if (i % 16 == 15)
            puts("");
    }
    fflush(stdout);
#endif
    int32_t r30 = unk6 << state->unk_shift;
    printf("r30: 0x%08X\n", r30); // 0x6400
    int32_t result[4][4];
    if (block_type == 1)
        r30 -= GetAot1(state, result, state->nest_data, state->h_nest_size, plane_idx);
    else
        r30 -= GetAotSum(state, result, block_type, state->nest_data, state->h_nest_size, plane_idx);
#if 1
    printf("r30: 0x%08X\n", r30); // 0x6400 - 0x3A00 => 0x2A00, not 0xFFFCC400
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
            printf("result[%u][%u]: 0x%08X\n", i, j, result[i][j]);
        }
#endif
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            uint32_t address = ((result[i][j] + r30) >> state->unk_shift);
            uint8_t byte = clipTable[address + 0x80];
            dst[i * stride + j] = byte;
        }
    }
}

static void printStackState(StackState *s)
{
    printf("{plane=%u, %p,%p,%p, %u,%u,%u,%u, %u}",
            s->plane_idx, s->ptr4, s->ptr8, s->ptrC,
            s->unk10, s->unk11, s->unk12, s->unk13, s->unk14);
}

static void IpicBlockDec(VideoState *state, void *present, uint32_t stride, StackState *stack_state)
{
    //printf(" IpicBlockDec(present=%p, stride=%u, ", present, stride);
    //printStackState(stack_state);
    //printf(")\n");
    if (stack_state->unk13 == 0)
    {
        uint8_t const *ptr = stack_state->ptr4;
        uint8_t r25 = ptr[1] & 0x77 ? stack_state->unk12 : ptr[0];
        ptr = stack_state->ptrC;
        uint8_t r24 = ptr[1] & 0x77 ? stack_state->unk12 : ptr[0];
        uint8_t r23 = stack_state->unk11 & 0x77 ? stack_state->unk12 : stack_state->unk10;
        WeightImBlock(present, stride, stack_state->unk12, r25, r24, stack_state->unk14, r23);
        stack_state->unk14 = stack_state->unk12;
    }
    else if (stack_state->unk13 == 8)
    {
        dcBlock(present, stride, stack_state->unk12);
        stack_state->unk14 = stack_state->unk12;
    }
    else
    {
        IntraAotBlock(state, present, stride, stack_state->unk12, stack_state->unk13, stack_state->plane_idx);
        stack_state->unk14 = stack_state->unk10;
    }
    stack_state->ptr4 += 2;
    stack_state->ptrC += 2;
}

static void IpicLineDec(VideoState *state, void *present, uint32_t stride, StackState *stack_state, uint16_t h_blocks)
{
    //printf("IpicLineDec(present=%p, stride=%u, ", present, stride);
    //printStackState(stack_state);
    //printf(", h_blocks=%u)\n", h_blocks);
    uint8_t *ptr = stack_state->ptr8;
    stack_state->unk10 = ptr[0];
    stack_state->unk11 = ptr[1];
    stack_state->unk14 = ptr[0];

    while (--h_blocks > 0)
    {
        stack_state->unk12 = stack_state->unk10;
        stack_state->unk13 = stack_state->unk11;
        stack_state->ptr8 += 2;
        ptr = stack_state->ptr8;
        stack_state->unk10 = ptr[0];
        stack_state->unk11 = ptr[1];
        IpicBlockDec(state, present, stride, stack_state);
        present += 4;
    }

    stack_state->unk12 = stack_state->unk10;
    stack_state->unk13 = stack_state->unk11;
    stack_state->ptr8 += 6;
    ptr = stack_state->ptr8;
    stack_state->unk10 = ptr[0];
    stack_state->unk11 = ptr[1];
    IpicBlockDec(state, present, stride, stack_state);

    stack_state->ptr4 += 4;
    stack_state->ptrC += 4;
}

static void IpicPlaneDec(VideoState *state, int plane_idx, void *present)
{
    printf("IpicPlaneDec(plane=%u, present=%p)\n", plane_idx, present);
    HVQPlaneDesc *plane = &state->planes[plane_idx];
    StackState stack_state;
    stack_state.plane_idx = plane_idx;
    stack_state.ptr4 = plane->ptr4;
    stack_state.ptr8 = plane->ptr4;
    stack_state.ptrC = plane->ptr4 + plane->h_blocks_safe * sizeof(uint16_t);
    int16_t v_blocks = plane->v_blocks;
    if (v_blocks > 0)
    {
        IpicLineDec(state, present, plane->width_in_samples, &stack_state, plane->h_blocks);
        present += plane->width_in_samples * 4;
        --v_blocks;
    }
    stack_state.ptr4 = plane->ptr4;
    while (v_blocks > 1)
    {
        IpicLineDec(state, present, plane->width_in_samples, &stack_state, plane->h_blocks);
        present += plane->width_in_samples * 4;
        --v_blocks;
    }
    if (v_blocks > 0)
    {
        stack_state.ptrC = stack_state.ptr8;
        IpicLineDec(state, present, plane->width_in_samples, &stack_state, plane->h_blocks);
    }
}

static void HVQM4DecodeIpic(SeqObj *seqobj, uint8_t const *frame, void *present)
{
    VideoState *state = seqobj->state;
    uint8_t scale = frame[0];
    state->unk_shift = frame[1];
    uint16_t foo = read16(frame + 4);
    uint16_t bar = read16(frame + 6);
    uint8_t const *data = frame + 0x48;
    frame += 8;
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
        setCode(&state->huffBuff[i].buf, data + read32(frame)); frame += 4;
    }
    readTree(&state->bufTree1[0], 0, 0);
    readTree(&state->bufTree2[0], 0, 0);
    readTree(&state->symbBuff[0], 1, scale);
    readTree(&state->bufTree0[0], 0, 2);

    state->boundB = +0x7F << scale;
    state->boundA = -0x80 << scale;

    Ipic_BasisNumDec(state);
    IpicDcvDec(state);
    MakeNest(state, foo, bar);

    for (int i = 0; i < 3; ++i)
    {
        IpicPlaneDec(state, i, present);
        present += state->planes[i].size_in_samples;
    }
}

static void HVQM4DecodeBpic(SeqObj *seqobj, uint8_t const *frame, void *present, void *past, void *future)
{
    // TODO
}

static void HVQM4DecodePpic(SeqObj *seqobj, uint8_t const *frame, void *present, void *past)
{
    HVQM4DecodeBpic(seqobj, frame, present, past, past);
}

static void decode_video(SeqObj *seqobj, FILE *infile, uint16_t frame_type, uint32_t frame_size)
{
    uint32_t pos = ftell(infile);
    uint8_t *frame = malloc(frame_size + 3 /* make ASan happy about decodeHuff() */);
    fread(frame, frame_size, 1, infile);
#if 0
    for (uint32_t i = 0; i < frame_size; ++i)
    {
        printf("%02X ", frame[i]);
        if (i % 16 == 15)
            puts("");
    }
    puts("");
#endif
    uint32_t frame_id = read32(frame);
    printf("display order within GOP: %u\n", frame_id);
    // HACK
    void *present = seqobj->state->yuvbuf, *past = present, *future = present;
    switch (frame_type)
    {
    case 0x10: puts("I frame"); HVQM4DecodeIpic(seqobj, frame + 4, present);               break;
    case 0x20: puts("P frame"); HVQM4DecodePpic(seqobj, frame + 4, present, past);         break;
    case 0x30: puts("B frame"); HVQM4DecodeBpic(seqobj, frame + 4, present, past, future); break;
    default:
        fprintf(stderr, "unknown video frame type 0x%x\n", frame_type);
        exit(EXIT_FAILURE);
    }
    fseek(infile, pos + frame_size, SEEK_SET);
    free(frame);
    puts("");

    if (frame_type == 0x10)
    {
        char name[50];
        static int i = 0; // HACK
        sprintf(name, "video_frame_%02u.ppm", i++);
        printf("got an I frame! writing it to %s now...\n", name);
        dumpYUV(seqobj, name);
    }
}

/* stream structure */

const char HVQM4_13_magic[16] = "HVQM4 1.3";
const char HVQM4_15_magic[16] = "HVQM4 1.5";

struct HVQM4_header
{
    enum
    {
        HVQM4_13,
        HVQM4_15,
    } version;

    uint32_t header_size;   /* 0x10-0x13 */
    uint32_t body_size;     /* 0x14-0x17 */
    uint32_t blocks;        /* 0x18-0x1B */
    uint32_t audio_frames;  /* 0x1C-0x1F */
    uint32_t video_frames;  /* 0x20-0x23 */
    uint32_t usec_per_frame;/* 0x24-0x27 (33366, 33367, 40000) */
    uint32_t duration;      /* 0x28-0x2B */
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
};

static void load_header(struct HVQM4_header *header, uint8_t *raw_header)
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
    header->audio_frames = read32(&raw_header[0x1C]);
    header->video_frames = read32(&raw_header[0x20]);
    header->usec_per_frame = read32(&raw_header[0x24]);
    header->duration = read32(&raw_header[0x28]);
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
    if (header->audio_frames != 0)
    {
        if (header->audio_srate == 0 || header->audio_frame_sz == 0 || header->audio_channels == 0)
        {
            fprintf(stderr, "expected nonzero audio srate and frame size\n");
            exit(EXIT_FAILURE);
        }
    }
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

void display_header(struct HVQM4_header *header)
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
    printf("Duration:    %"PRIu32" (?)\n", header->duration);
    printf("Resolution:  %"PRIu32" x %"PRIu32"\n", header->hres, header->vres);
    printf("µs/frame:    %"PRIu32"\n", header->usec_per_frame);
    printf("%d Blocks\n", header->blocks);
    printf("%d Video frames\n", header->video_frames);
    if (header->audio_frames)
    {
        printf("%d Audio frames\n", header->audio_frames);
        printf("Sample rate: %"PRIu32" Hz\n", header->audio_srate);
        printf("Audio frame size: 0x%"PRIx32"\n", header->audio_frame_sz);
        printf("Audio channels: %u\n", header->audio_channels);
    } else {
        printf("No audio!\n");
    }
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

int main(int argc, char **argv)
{
    printf("h4m 'HVQM4 1.3/1.5' audio decoder 0.3 by hcs\n\n");
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
    struct HVQM4_header header;
    load_header(&header, raw_header);
    display_header(&header);

    if (header.audio_frames == 0)
    {
        fprintf(stderr, "this video contains no audio!\n");
        exit(EXIT_FAILURE);
    }

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

    SeqObj seqobj;
    HVQM4InitDecoder();
    HVQM4InitSeqObj(&seqobj, (VideoInfo*)&header.hres);
    VideoState *state = malloc(HVQM4BuffSize(&seqobj));
    HVQM4SetBuffer(&seqobj, state);

    /* parse blocks */
    uint32_t block_count = 0;
    uint32_t total_aud_frames = 0;
    uint32_t total_vid_frames = 0;
    uint32_t last_block_start = ftell(infile);
    uint32_t total_sample_count = 0;
    while (block_count < header.blocks)
    {
        const long block_start = ftell(infile);
        expect32(ftell(infile) - last_block_start, infile);
        last_block_start = block_start;
        const uint32_t expected_block_size = get32(infile);
        const uint32_t expected_aud_frame_count = get32(infile);
        const uint32_t expected_vid_frame_count = get32(infile);
        expect32(0x01000000, infile);   /* EOS marker? */
        const long data_start = ftell(infile);

        block_count ++;
#ifdef VERBOSE_PRINT
        printf("\n\nblock %d starts at 0x%lx, length 0x%"PRIx32"\n", (int)block_count, block_start, expected_block_size);
#endif

        /* parse frames */
        struct audio_state audio_state;
        int first_vid=1, first_aud=1;
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
            printf("\nframe id 0x%"PRIx16",0x%"PRIx16" ",frame_id1,frame_id2);
            printf("size 0x%"PRIx32"\n", frame_size);
#endif

            if (frame_id1 == 1)
            {
                /* video */
                first_vid = 0;
                vid_frame_count ++;
                total_vid_frames ++;
#ifdef VERBOSE_PRINT
                printf("video frame %d/%d (%d)\n", (int)vid_frame_count, (int)expected_vid_frame_count, (int)total_vid_frames);
#endif
                decode_video(&seqobj, infile, frame_id2, frame_size);
            }
            else if (frame_id1 == 0 &&
                ((first_aud && ( frame_id2 == 3 || frame_id2 == 1)) ||
                (!first_aud && frame_id2 == 2)))

            {
                /* audio */
                const long audio_started = ftell(infile);
                const uint32_t samples = get32(infile);
                decode_audio(&audio_state, first_aud, samples, infile, outfile, header.audio_channels);
                block_sample_count += samples;
                aud_frame_count ++;
                total_aud_frames ++;
#ifdef VERBOSE_PRINT
                printf("0x%lx: audio frame %d/%d (%d) (%d samples)\n", (unsigned long)audio_started, (int)aud_frame_count, (int)expected_aud_frame_count, (int)total_aud_frames, samples);
#endif
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
            }
            else
            {
                fprintf(stderr, "unexpected frame id at %08lx\n", (unsigned long)ftell(infile));
                exit(EXIT_FAILURE);
                printf("unexpected frame id %d %d at %08lx\n", frame_id1, frame_id2, (unsigned long)(ftell(infile)-8));
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
