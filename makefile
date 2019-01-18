#SAMPLE = samples/LOGOS.h4m
SAMPLE = samples/pikminS.h4m
REFERENCE = reference/

build_emu:
	LD_LIBRARY_PATH=toolchain/lib toolchain/bin/powerpc-linux-gcc -Wall -Wextra -Wno-unused-function -g -fno-omit-frame-pointer -static h4m_audio_decode.c -o h4m_audio_decode
	rm -f output/*.ppm

emu: build_emu
	qemu-ppc h4m_audio_decode $(SAMPLE) foo.wav

debug: build_emu
	qemu-ppc -g 1234 h4m_audio_decode $(SAMPLE) foo.wav &
	toolchain/bin/powerpc-linux-gdb -ex 'target remote localhost:1234' -ex c h4m_audio_decode

native:
	clang -m32 -march=native -O2 -funroll-loops -Wall -Wextra h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f output/*.ppm
	time -p ./h4m_audio_decode $(SAMPLE) foo.wav
	diff -q output/ $(REFERENCE)

native64:
	clang -march=native -O2 -funroll-loops -Wall -Wextra h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f output/*.ppm
	time -p ./h4m_audio_decode $(SAMPLE) foo.wav
	diff -q output/ $(REFERENCE)


clean:
	rm -f h4m_audio_decode output/*.ppm

.PHONY: all clean
