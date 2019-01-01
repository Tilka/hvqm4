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
    void *f;
    size_t f_size;
    Elf32_Shdr *symtab;
    Elf32_Shdr *strtab;
    Elf32_Shdr *relatext;
    void *trampolines;
} yolo_lib;

void yolo_open(yolo_lib *lib, const char *path)
{
    //printf("loading %s...\n", path);
    FILE *f = fopen(path, "rb");
    if (!f)
    {
        perror("fopen");
        exit(1);
    }

    fseek(f, 0, SEEK_END);
    lib->f_size = ftell(f);
    fseek(f, 0, SEEK_SET);

    lib->f = mmap(NULL, lib->f_size, PROT_READ, MAP_PRIVATE, fileno(f), 0);
    if (lib->f == MAP_FAILED)
    {
        perror("mmap file");
        exit(1);
    }
    fclose(f);

    void *max_addr = 0;
    Elf32_Ehdr *ehdr = lib->f;
    long page_size = sysconf(_SC_PAGESIZE);
    for (size_t i = 0; i < ehdr->e_phnum; ++i)
    {
        Elf32_Phdr *phdr = lib->f + ehdr->e_phoff + i * sizeof(Elf32_Phdr);
        if (phdr->p_type != PT_LOAD)
            continue;

        void *addr = (void*)phdr->p_vaddr;
        void *aligned_addr = (void*)((long)addr & ~(page_size - 1));
        size_t size = phdr->p_memsz;
        size_t aligned_size = size + (addr - aligned_addr);
        int prot = 0;
        if (phdr->p_flags & PF_R)
            prot |= PROT_READ;
        if (phdr->p_flags & PF_W)
            prot |= PROT_WRITE;
        if (phdr->p_flags & PF_X)
            prot |= PROT_EXEC;
#if 0
        printf("mmap: addr 0x%08X aaddr 0x%08X size 0x%08X asize 0x%08X fsize 0x%08X offset 0x%08X %c%c%c\n", addr, aligned_addr, size, aligned_size, phdr->p_filesz, phdr->p_offset,
               "r-"[!(phdr->p_flags & PF_R)],
               "w-"[!(phdr->p_flags & PF_W)],
               "x-"[!(phdr->p_flags & PF_X)]);
#endif
        // mmap requires page-aligned file offsets, so we can't just map the file :(
        if (mmap(aligned_addr, aligned_size, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_ANONYMOUS | MAP_FIXED | MAP_PRIVATE, -1, 0) != aligned_addr)
        {
            perror("mmap");
            exit(1);
        }
        max_addr = addr + size > max_addr ? addr + size : max_addr;
    }

    // map one adjacent page to store patch trampolines
    lib->trampolines = (void*)((long)max_addr & ~(page_size - 1));
    if (mmap(lib->trampolines, page_size, PROT_READ | PROT_WRITE | PROT_EXEC, MAP_ANONYMOUS | MAP_FIXED | MAP_PRIVATE, -1, 0) != lib->trampolines)
    {
        perror("mmap trampolines");
        exit(1);
    }

    // need to first mmap all segments (in case we have overlaps due to alignment) and *then* copy the contents,
    // otherwise the kernel will zero-out the overlap
    for (size_t i = 0; i < ehdr->e_phnum; ++i)
    {
        Elf32_Phdr *phdr = lib->f + ehdr->e_phoff + i * sizeof(Elf32_Phdr);
        if (phdr->p_type != PT_LOAD)
            continue;

        void *addr = (void*)phdr->p_vaddr;
        size_t size = phdr->p_memsz;
        memcpy(addr, lib->f + phdr->p_offset, phdr->p_filesz);
#if 0
        long page_size = sysconf(_SC_PAGESIZE);
        void *aligned_addr = (void*)((long)addr & ~(page_size - 1));
        size_t aligned_size = size + (addr - aligned_addr);
        int prot = 0;
        if (phdr->p_flags & PF_R)
            prot |= PROT_READ;
        if (phdr->p_flags & PF_W)
            prot |= PROT_WRITE;
        if (phdr->p_flags & PF_X)
            prot |= PROT_EXEC;
        if (mprotect(aligned_addr, aligned_size, prot) != 0)
        {
            perror("mprotect");
            exit(1);
        }
#endif
    }

    Elf32_Shdr *shstr = lib->f + ehdr->e_shoff + ehdr->e_shstrndx * sizeof(Elf32_Shdr);
    for (size_t i = 0; i < ehdr->e_shnum; ++i)
    {
        Elf32_Shdr *shdr = lib->f + ehdr->e_shoff + i * sizeof(Elf32_Shdr);
        if (shdr->sh_type == SHT_SYMTAB)
            lib->symtab = shdr;
        if (shdr->sh_type == SHT_STRTAB)
            lib->strtab = shdr;
        const char *shname = lib->f + shstr->sh_offset + shdr->sh_name;
        if (strcmp(shname, ".rela.text") == 0)
            lib->relatext = shdr;
    }
}

size_t yolo_sym_index(yolo_lib *lib, const char *symbol)
{
    size_t sym_count = lib->symtab->sh_size / sizeof(Elf32_Sym);
    for (size_t i = 0; i < sym_count; ++i)
    {
        Elf32_Sym *sym = lib->f + lib->symtab->sh_offset + i * sizeof(Elf32_Sym);
        const char *name = lib->f + lib->strtab->sh_offset + sym->st_name;
        if (strcmp(name, symbol) == 0)
            return i;
    }
    fprintf(stderr, "couldn't find symbol '%s'\n", symbol);
    exit(EXIT_FAILURE);
}

void *yolo_sym(yolo_lib *lib, const char *symbol)
{
    size_t i = yolo_sym_index(lib, symbol);
    Elf32_Sym *sym = lib->f + lib->symtab->sh_offset + i * sizeof(Elf32_Sym);
    return sym->st_value;
}

// patch all references to symbol so that they point to target
static void yolo_patch(yolo_lib *lib, const char *symbol, void *target)
{
    Elf32_Ehdr *ehdr = lib->f + lib->relatext->sh_offset;
    size_t sym_index = yolo_sym_index(lib, symbol);
    size_t rel_count = lib->relatext->sh_size / sizeof(Elf32_Rela);
    for (size_t i = 0; i < rel_count; ++i)
    {
        Elf32_Rela *rel = lib->f + lib->relatext->sh_offset + i * sizeof(Elf32_Rela);
        if (ELF32_R_SYM(rel->r_info) == sym_index)
        {
            printf("need to patch reference to %s at 0x%08X\n", symbol, rel->r_offset);
            uint32_t *ptr = (uint32_t*)rel->r_offset;
            size_t diff = lib->trampolines - rel->r_offset;
            if (diff > 0x3FFFFFF)
            {
                show_mappings();
                fputs("R_PPC_REL24 target out-of-range\n", stderr);
                exit(1);
            }
            uint32_t mask = 0x3FFFFFC;
            *ptr = (*ptr & ~mask) | (diff & mask);
            uint32_t *code = lib->trampolines;
            // lis r0, hi
            code[0] = 0x3C000000 | ((long)target >> 16);
            // ori r0, lo
            code[1] = 0x60000000 | ((long)target & 0xFFFF);
            // mtctr
            code[2] = 0x7C0903A6;
            // bctr
            code[3] = 0x4E800420;
            lib->trampolines += 4 * sizeof(uint32_t);
        }
    }
}

static void yolo_close(yolo_lib *lib)
{
    munmap(lib->f, lib->f_size);
}

static void load_library(const char *lib_path)
{
    yolo_lib lib;
    yolo_open(&lib, lib_path);

    // locate symbols
#define SYMBOLT(x, y) p##x = yolo_sym(&lib, #x);
#include "symbols.inc"
#undef SYMBOLT

    // patch function calls
#define PATCH(x) yolo_patch(&lib, #x, x);
    PATCH(GetMCAot1)
    PATCH(GetMCAotSum)
    PATCH(_MotionComp)
    PATCH(decodeSOvfSym)
    PATCH(decodeHuff)
#undef PATCH

    yolo_close(&lib);
}
