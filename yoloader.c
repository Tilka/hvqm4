#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <elf.h>
#include <string.h>
#include <unistd.h>

static void show_mappings()
{
    puts("mappings:");
    FILE *f = fopen("/proc/self/maps", "r");
    int c;
    while((c = fgetc(f)) != EOF)
        printf("%c", c);
    puts("");
}

typedef struct
{
    FILE *f;
    Elf32_Ehdr ehdr;
} yolo_lib;

void yolo_open(yolo_lib *lib, const char *path)
{
    //printf("loading %s...\n", path);
    FILE *f = lib->f = fopen(path, "rb");
    if (!f)
    {
        perror("fopen");
        exit(1);
    }

    Elf32_Ehdr *ehdr = &lib->ehdr;
    if (fread(ehdr, sizeof(Elf32_Ehdr), 1, f) != 1)
    {
        perror("fread ehdr");
        exit(1);
    }
    if (fseek(f, ehdr->e_phoff, SEEK_SET) != 0)
    {
        perror("fseek phoff");
        exit(1);
    }

    for (unsigned i = 0; i < ehdr->e_phnum; ++i)
    {
        Elf32_Phdr phdr;
        if (fread(&phdr, sizeof(Elf32_Phdr), 1, f) != 1)
        {
            perror("fread phdr");
            exit(1);
        }
        if (phdr.p_type != PT_LOAD)
            continue;
        long next_phdr = ftell(f);

        long page_size = sysconf(_SC_PAGESIZE);
        void *addr = (void*)phdr.p_vaddr;
        void *aligned_addr = (void*)((long)addr & ~(page_size - 1));
        size_t size = phdr.p_memsz;
        size_t aligned_size = size + (addr - aligned_addr);
        int prot = 0;
        if (phdr.p_flags & PF_R)
            prot |= PROT_READ;
        if (phdr.p_flags & PF_W)
            prot |= PROT_WRITE;
        if (phdr.p_flags & PF_X)
            prot |= PROT_EXEC;
#if 0
        printf("mmap: addr 0x%08X aaddr 0x%08X size 0x%08X asize 0x%08X fsize 0x%08X offset 0x%08X %c%c%c\n", addr, aligned_addr, size, aligned_size, phdr.p_filesz, phdr.p_offset,
               "r-"[!(phdr.p_flags & PF_R)],
               "w-"[!(phdr.p_flags & PF_W)],
               "x-"[!(phdr.p_flags & PF_X)]);
#endif
        // mmap requires page-aligned file offsets, so we can't just map the file :(
        if (mmap(aligned_addr, aligned_size, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_FIXED | MAP_PRIVATE, -1, 0) != aligned_addr)
        {
            perror("mmap");
            exit(1);
        }
        if (phdr.p_filesz)
        {
            if (fseek(f, phdr.p_offset, SEEK_SET) != 0)
            {
                perror("fseek content");
                exit(1);
            }
            size_t res = fread(addr, phdr.p_filesz, 1, f);
            if (res != 1)
            {
                fprintf(stderr, "res: %zu eof: %u error: %u\n", res, feof(f), ferror(f));
                perror("fread content");
                exit(1);
            }
            if (fseek(f, next_phdr, SEEK_SET) != 0)
            {
                perror("fseek next phdr");
                exit(1);
            }
        }
        if (mprotect(aligned_addr, aligned_size, prot) != 0)
        {
            perror("mprotect");
            exit(1);
        }
    }
}

void *yolo_sym(yolo_lib *lib, const char *symbol)
{
#if 0
    FILE *f = lib->f;
    if (fseek(f, lib->ehdr.e_shoff, SEEK_SET) != 0)
    {
        perror("fseek shoff");
        exit(1);
    }
#endif

    // TODO lol
    if (strcmp(symbol, "HVQM4InitDecoder") == 0)
        return (void*)0x801B5DF0;
    if (strcmp(symbol, "HVQM4InitSeqObj") == 0)
        return (void*)0x801B5E10;
    if (strcmp(symbol, "HVQM4BuffSize") == 0)
        return (void*)0x801B5E34;
    if (strcmp(symbol, "HVQM4SetBuffer") == 0)
        return (void*)0x801B5EB8;
    if (strcmp(symbol, "HVQM4DecodeIpic") == 0)
        return (void*)0x801B6144;
    if (strcmp(symbol, "HVQM4DecodePpic") == 0)
        return (void*)0x801B6340;
    if (strcmp(symbol, "HVQM4DecodeBpic") == 0)
        return (void*)0x801B638C;
    if (strcmp(symbol, "getMVector") == 0)
        return (void*)0x801B5CF4;
    if (strcmp(symbol, "MCBlockDecDCNest") == 0)
        return (void*)0x801B54D8;
    if (strcmp(symbol, "MCBlockDecMCNest") == 0)
        return (void*)0x801B5330;
    if (strcmp(symbol, "PrediAotBlock") == 0)
        return (void*)0x801B45A8;
    if (strcmp(symbol, "IpicPlaneDec") == 0)
        return (void*)0x801B38C0;
    if (strcmp(symbol, "IpicLineDec") == 0)
        return (void*)0x801B37E4;
    if (strcmp(symbol, "IpicBlockDec") == 0)
        return (void*)0x801B3680;
    if (strcmp(symbol, "initMCHandler") == 0)
        return (void*)0x801B39E0;
    if (strcmp(symbol, "spread_PB_descMap") == 0)
        return (void*)0x801B5854;
    if (strcmp(symbol, "resetMCHandler") == 0)
        return (void*)0x801B3AE0;
    if (strcmp(symbol, "MotionComp") == 0)
        return (void*)0x801B3B38;
    fputs("TODO\n", stderr);
    exit(1);
}

static void bla()
{
    fputs("called an uninitialized function pointer\n", stderr);
    exit(1);
}

void (*pHVQM4InitDecoder)() = bla;
void (*pHVQM4InitSeqObj)() = bla;
size_t (*pHVQM4BuffSize)() = bla;
void (*pHVQM4SetBuffer)() = bla;
void (*pHVQM4DecodeIpic)() = bla;
void (*pHVQM4DecodePpic)() = bla;
void (*pHVQM4DecodeBpic)() = bla;
void (*pgetMVector)() = bla;
void (*pMCBlockDecDCNest)() = bla;
void (*pMCBlockDecMCNest)() = bla;
void (*pPrediAotBlock)() = bla;
void (*pIpicPlaneDec)() = bla;
void (*pIpicLineDec)() = bla;
void (*pIpicBlockDec)() = bla;
void (*pinitMCHandler)() = bla;
void (*pspread_PB_descMap)() = bla;
void (*presetMCHandler)() = bla;
void (*pMotionComp)() = bla;

static void load_library(const char *lib_path)
{
    yolo_lib rm_dll;
    yolo_open(&rm_dll, lib_path);
    pHVQM4InitDecoder = yolo_sym(&rm_dll, "HVQM4InitDecoder");
    pHVQM4InitSeqObj  = yolo_sym(&rm_dll, "HVQM4InitSeqObj");
    pHVQM4BuffSize    = yolo_sym(&rm_dll, "HVQM4BuffSize");
    pHVQM4SetBuffer   = yolo_sym(&rm_dll, "HVQM4SetBuffer");
    pHVQM4DecodeIpic  = yolo_sym(&rm_dll, "HVQM4DecodeIpic");
    pHVQM4DecodePpic  = yolo_sym(&rm_dll, "HVQM4DecodePpic");
    pHVQM4DecodeBpic  = yolo_sym(&rm_dll, "HVQM4DecodeBpic");

    pgetMVector       = yolo_sym(&rm_dll, "getMVector");
    pMCBlockDecDCNest = yolo_sym(&rm_dll, "MCBlockDecDCNest");
    pMCBlockDecMCNest = yolo_sym(&rm_dll, "MCBlockDecMCNest");
    pPrediAotBlock    = yolo_sym(&rm_dll, "PrediAotBlock");
    pIpicPlaneDec     = yolo_sym(&rm_dll, "IpicPlaneDec");
    pIpicLineDec      = yolo_sym(&rm_dll, "IpicLineDec");
    pIpicBlockDec     = yolo_sym(&rm_dll, "IpicBlockDec");
    pinitMCHandler    = yolo_sym(&rm_dll, "initMCHandler");
    pspread_PB_descMap = yolo_sym(&rm_dll, "spread_PB_descMap");
    presetMCHandler   = yolo_sym(&rm_dll, "resetMCHandler");
    pMotionComp       = yolo_sym(&rm_dll, "MotionComp");
}

#ifndef YOLO_INCLUDE
int main(int argc, char **argv)
{
    load_library("RM_DLL.elf");

    pHVQM4InitDecoder();
    SeqObj seqobj;
    char header[0x44];
    FILE *h4m = fopen(argv[1], "rb");
    fread(&header, sizeof(header), 1, h4m);
    pHVQM4InitSeqObj(&seqobj, &header[0x34]);
    size_t bufsize = pHVQM4BuffSize(seqobj);
    void *workbuf = malloc(bufsize);
    pHVQM4SetBuffer(&seqobj, workbuf);



    puts("success");
}
#endif
