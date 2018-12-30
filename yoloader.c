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
    char const **section_names;
    Elf32_Shdr *symtab;
    Elf32_Shdr *strtab;
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

    Elf32_Ehdr *ehdr = lib->f;
    for (unsigned i = 0; i < ehdr->e_phnum; ++i)
    {
        Elf32_Phdr *phdr = lib->f + ehdr->e_phoff + i * sizeof(Elf32_Phdr);
        if (phdr->p_type != PT_LOAD)
            continue;

        long page_size = sysconf(_SC_PAGESIZE);
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
        memcpy(addr, lib->f + phdr->p_offset, phdr->p_filesz);
        if (mprotect(aligned_addr, aligned_size, prot) != 0)
        {
            perror("mprotect");
            exit(1);
        }
    }

    for (unsigned i = 0; i < ehdr->e_shnum; ++i)
    {
        Elf32_Shdr *shdr = lib->f + ehdr->e_shoff + i * sizeof(Elf32_Shdr);
        if (shdr->sh_type == SHT_SYMTAB)
            lib->symtab = shdr;
        if (shdr->sh_type == SHT_STRTAB)
            lib->strtab = shdr;
    }
}

void *yolo_sym(yolo_lib *lib, const char *symbol)
{
    char const *str_start = lib->f + lib->strtab->sh_offset;
    char const *str = str_start;
    while (str - str_start < lib->strtab->sh_size)
    {
        if (strcmp(str, symbol) == 0)
        {
            for (unsigned j = 0; j < lib->symtab->sh_size / sizeof(Elf32_Sym); ++j)
            {
                Elf32_Sym *sym = lib->f + lib->symtab->sh_offset + j * sizeof(Elf32_Sym);
                if (sym->st_name == str - str_start)
                    return sym->st_value;
            }
            fprintf(stderr, "couldn't find symbol '%s'\n", symbol);
            exit(EXIT_FAILURE);
        }
        str += strlen(str) + 1;
    }
    fprintf(stderr, "couldn't find symbol '%s'\n", symbol);
    exit(EXIT_FAILURE);
}

static void yolo_close(yolo_lib *lib)
{
    munmap(lib->f, lib->f_size);
}

static void bla()
{
    fputs("called an uninitialized function pointer\n", stderr);
    exit(1);
}

#define SYMBOLT(x, T) T (*p##x)() = bla;
#include "symbols.inc"
#undef SYMBOLT

static void load_library(const char *lib_path)
{
    yolo_lib rm_dll;
    yolo_open(&rm_dll, lib_path);
#define SYMBOLT(x, y) p##x = yolo_sym(&rm_dll, #x);
#include "symbols.inc"
#undef SYMBOLT
}
