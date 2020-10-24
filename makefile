#MOVIE = LOGOS.h4m
#MOVIE = PKLOGO.h4m
#MOVIE = pikminS.h4m
#MOVIE = theix.h4m
MOVIE = Treyarch_Logo.h4m
#MOVIE = crap starfox.h4m
SAMPLE = "samples/$(MOVIE)"
REFERENCE = "reference/$(MOVIE)"
DIFF = diff -q --binary --speed-large-files

build_emu:
	LD_LIBRARY_PATH=toolchain/lib toolchain/bin/powerpc-linux-gcc -Og -Wall -Wextra -Wno-unused-function -g -fno-omit-frame-pointer -static h4m_audio_decode.c -o h4m_audio_decode
	rm -f output/*.ppm

emu: build_emu
	qemu-ppc h4m_audio_decode $(SAMPLE) foo.wav
	$(DIFF) output/ $(REFERENCE)

debug: build_emu
	qemu-ppc -g 1234 h4m_audio_decode $(SAMPLE) foo.wav &
	toolchain/bin/powerpc-linux-gdb -ex 'target remote localhost:1234' -ex c h4m_audio_decode

native:
	clang -m32 -march=native -O2 -funroll-loops -Wall -Wextra h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f output/*.ppm
	time -p ./h4m_audio_decode $(SAMPLE) foo.wav
	$(DIFF) -q output/ $(REFERENCE)

native64:
	clang -march=native -O2 -funroll-loops -Wall -Wextra h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f output/*.ppm
	time -p ./h4m_audio_decode $(SAMPLE) foo.wav
	$(DIFF) -q output/ $(REFERENCE)


clean:
	rm -f h4m_audio_decode output/*.ppm

.PHONY: all clean
