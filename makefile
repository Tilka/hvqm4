build_emu:
	toolchain/bin/powerpc-linux-gcc -Wall -Wextra -Wno-unused-function -g -fno-omit-frame-pointer -static h4m_audio_decode.c -o h4m_audio_decode
	rm -f output/*.ppm

emu: build_emu
	qemu-ppc h4m_audio_decode samples/LOGOS.h4m foo.wav

debug: build_emu
	qemu-ppc -g 1234 h4m_audio_decode samples/LOGOS.h4m foo.wav &
	toolchain/bin/powerpc-linux-gdb -ex 'target remote localhost:1234' -ex c h4m_audio_decode

native:
	clang -Wall -Wextra -Wno-unused-function -g -fsanitize=address -fno-omit-frame-pointer -m32 h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f output/*.ppm
	./h4m_audio_decode samples/LOGOS.h4m foo.wav

clean:
	rm -f h4m_audio_decode output/*.ppm

.PHONY: all clean
